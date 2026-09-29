#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double sum(double *x, int dim) {
	double s;
    int i;
    s = 0.0;
    for (i=0; i<dim; ++i) s+=x[i];
    	return(s);
}


void doDL(double *x, double *le, double *p1, double *p2, int *dim_le, 
                                double *dl_values){
	double m1, m2, k2;
    double d1, d2, d3, d4, k1, z;
    int i;
        
    d1 = x[0];
    d2 = x[1];
    d3 = x[2];
    d4 = x[3];
    k1 = x[4];
    z = x[5];
    
    m1 = d1 + 0.5*d2;
    m2 = d1 + d2 + d3 + 0.5*d4;
    k2 = z - k1;

    for (i=0; i< (*dim_le); i++){
		dl_values[i] = k1/(1+exp(-log(pow(*p1,2.0))*(le[i]-m1)/d2)) + k2/(1+exp(-log(pow(*p2,2.0))*(le[i]-m2)/d4));
    }
}

double rnormtrunc(double mu, double sigma, double low, double high){
	double temp;
	int maxit, i;
	
	temp = -999;
  	maxit = 1000;
  	i = 0;
  	GetRNGstate();
  	while((temp<low || temp>high) && i <= maxit) {
    	temp = rnorm(mu, sigma);
     	i++;
  	}
  	if (i > maxit) {
  		if(temp<low) temp = low;
  		else temp = high;
  	}
  	PutRNGstate();
  	return(temp);
}


void dnormtrunc(double *x, double *mu, double *sigma, 
		double low, double high, int dim_out, double *out){
	int i;
	for (i=0; i< dim_out; i++) {
		if(x[i] < low || x[i] > high) out[i] = 0;
		else 
  		out[i] = dnorm(x[i],mu[i],sigma[i], 0)/(pnorm(high,mu[i],sigma[i],1,0)-pnorm(low,mu[i],sigma[i],1,0));
  	}
	return;
}


void dologdensityTrianglekz(double *x, double *mu, double *sigma, 
			double *low, double *up, int *par_idx, double *dlpars, double *p1, double *p2,
			double *le, int *lidx, double *dct, double *loess_sd, double *logdens) {
	double dl[*lidx], dens[*lidx], dnt[1];
	double s;
	int i, param_index;
	param_index = *par_idx;
	
	if (param_index < 1) error("Wrong parameter index: %i", param_index);
	dlpars[param_index-1] = *x;
	doDL(dlpars, le, p1, p2, lidx, dl);
	s = 0;
	for (i=0; i< (*lidx); i++){
		dens[i] = dnorm(dct[i], dl[i], loess_sd[i], 0);
		/*Rprintf("\n%f, %f %f %f", dens[i], dct[i], dl[i], loess_sd[i]);*/
		if(dens[i] < 1e-100) dens[i] = 1e-100;
		s = s+ log(dens[i]);
	}
	dnormtrunc(x, mu, sigma, *low, *up, 1, dnt);
	logdens[0] = s + log(dnt[0]);
	/*Rprintf("\ns=%f dnt=%f res = %f", s, dnt[0], logdens[0]);*/
	return;
}

/* Log-density of one country-specific DL parameter (index idx, 1-based),
   evaluated at x with all other parameters taken from dlpars (work array). */
static double logdens_country_par(double x, int idx, double mean, double sd, double low, double up,
			double *dlpars, double *p1, double *p2, double *le, int n, double *dct, double *sdv) {
	double logdens;
	dologdensityTrianglekz(&x, &mean, &sd, &low, &up, &idx, dlpars, p1, p2, le, &n, dct, sdv, &logdens);
	return(logdens);
}

/* C version of the R function slice.sampling() applied to logdens_country_par.
   Draws random numbers in the same order as the R version. */
static double slice_sample_country_par(double x0, int idx, double width, double low, double up, 
			double mean, double sd, double *dlpars, double *p1, double *p2, 
			double *le, int n, double *dct, double *sdv) {
	int maxit = 50, i;
	double z, L, R, J, K, x1;
	
	z = logdens_country_par(x0, idx, mean, sd, low, up, dlpars, p1, p2, le, n, dct, sdv) - rexp(1.0);
	L = x0 - runif(0, width);
	R = L + width;
	J = floor(runif(0, maxit));
	K = (maxit-1) - J;
	while (J > 0 && L > low && 
			logdens_country_par(L, idx, mean, sd, low, up, dlpars, p1, p2, le, n, dct, sdv) > z) {
		L = L - width;
		J = J - 1;
	}
	while (K > 0 && R < up && 
			logdens_country_par(R, idx, mean, sd, low, up, dlpars, p1, p2, le, n, dct, sdv) > z) {
		R = R + width;
		K = K - 1;
	}
	if (L < low) L = low;
	if (R > up) R = up;
	if (L > R) return(x0);
	for (i = 1; i <= maxit; i++) {
		x1 = runif(L, R);
		if (z <= logdens_country_par(x1, idx, mean, sd, low, up, dlpars, p1, p2, le, n, dct, sdv))
			return(x1);
		if (x1 < x0) L = x1;
		else R = x1;
	}
	error("Problem in slice sampling");
	return(x0);
}

/* sum of dlx[0..3] without element i; long double to match R's sum() */
static double sum_Triangle_without(double *dlx, int i) {
	long double s = 0.0;
	int j;
	for (j = 0; j < 4; j++) if (j != i) s += dlx[j];
	return((double) s);
}

/* Update of Triangle.c (4 pars), k.c and z.c for one country via slice sampling. 
   All par vectors have length 6 in the order Triangle.c[1:4], k.c, z.c.
   sdv is omega * loess SD. Returns the updated 6 values. */
SEXP doTrianglekzcUpdate(SEXP scur, SEXP smean, SEXP ssd, SEXP slow, SEXP sup, SEXP swidth, 
			SEXP ssumlim, SEXP sp1, SEXP sp2, SEXP sle, SEXP sdct, SEXP ssdv) {
	double *cur = REAL(scur), *mean = REAL(smean), *sd = REAL(ssd), *low = REAL(slow),
			*up = REAL(sup), *width = REAL(swidth), *sumlim = REAL(ssumlim),
			*p1 = REAL(sp1), *p2 = REAL(sp2), *le = REAL(sle), *dct = REAL(sdct), *sdv = REAL(ssdv);
	int n = LENGTH(sle), i, ntries;
	double dlx[6], work[6], prop[4], lo, hi, s, sT;
	long double sTl;
	SEXP res;
	
	if (LENGTH(scur) != 6 || LENGTH(smean) != 6 || LENGTH(ssd) != 6 || LENGTH(slow) != 6 ||
			LENGTH(sup) != 6 || LENGTH(swidth) != 6 || LENGTH(ssumlim) != 2)
		error("Wrong length of parameter vectors.");
	if (LENGTH(sdct) != n || LENGTH(ssdv) != n) error("Inconsistent data lengths.");
	
	GetRNGstate();
	for (i = 0; i < 6; i++) dlx[i] = cur[i];
	for (i = 0; i < 4; i++) prop[i] = 0;
	ntries = 1;
	while (ntries <= 50) {
		for (i = 0; i < 4; i++) {
			s = sum_Triangle_without(dlx, i);
			lo = fmin2(fmax2(low[i], sumlim[0] - s), cur[i]);
			hi = fmax2(fmin2(up[i], sumlim[1] - s), cur[i]);
			memcpy(work, dlx, 6*sizeof(double));
			prop[i] = slice_sample_country_par(cur[i], i+1, width[i], lo, hi, mean[i], sd[i], 
								work, p1, p2, le, n, dct, sdv);
			dlx[i] = prop[i];
		}
		sTl = 0.0;
		for (i = 0; i < 4; i++) sTl += prop[i];
		sT = (double) sTl;
		if (sT <= sumlim[1] && sT >= sumlim[0]) break;
		for (i = 0; i < 6; i++) dlx[i] = cur[i];
		ntries++;
	}
	PROTECT(res = allocVector(REALSXP, 6));
	for (i = 0; i < 4; i++) REAL(res)[i] = prop[i];
	/* k.c and z.c */
	for (i = 4; i < 6; i++) {
		memcpy(work, dlx, 6*sizeof(double));
		dlx[i] = slice_sample_country_par(cur[i], i+1, width[i], low[i], up[i], mean[i], sd[i],
								work, p1, p2, le, n, dct, sdv);
		REAL(res)[i] = dlx[i];
	}
	PutRNGstate();
	UNPROTECT(1);
	return(res);
}

/* Double-logistic function for n parameter sets (rows of the n x 6 matrix sx),
   each evaluated at its own e0 value (sle, length n). */
SEXP doDLmulti(SEXP sx, SEXP sle, SEXP sp1, SEXP sp2) {
	int n = LENGTH(sle), i, j, one = 1;
	double *x = REAL(sx), *le = REAL(sle), pars[6];
	SEXP res;
	
	if (LENGTH(sx) != 6*n) error("Parameter matrix must have 6 columns and as many rows as there are e0 values.");
	PROTECT(res = allocVector(REALSXP, n));
	for (i = 0; i < n; i++) {
		for (j = 0; j < 6; j++) pars[j] = x[i + j*n];
		doDL(pars, &le[i], REAL(sp1), REAL(sp2), &one, &REAL(res)[i]);
	}
	UNPROTECT(1);
	return(res);
}
