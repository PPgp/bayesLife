#include <R.h>
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

void get_z_gibbs_pars(int *ncountries, int *ntime, double *e0matrix, double *loess, 
						double *Triangle_c, double *k_c, double *p1, double *p2,
						double *d) {
	int i, j, nrow, ncol;
	double A, B, y, x, m1, m2, lp1, lp2, b, c, h, oneplusexpy, oneplusexpx, tmp, a1;
	ncol = *ncountries;
	a1 = 0;
	b = 0;
	/*c = 0;*/
	lp1 = log(pow(*p1,2.0));
	lp2 = log(pow(*p2,2.0));
	for (j=0; j<ncol; j++) {
		nrow = ntime[j];
		m2 = Triangle_c[j*4] + Triangle_c[1 + j*4] + Triangle_c[2 + j*4] + 0.5*Triangle_c[3 + j*4];
		m1 = Triangle_c[j*4] + 0.5*Triangle_c[1 + j*4];
		for (i=0; i<(nrow-1); i++) {
			h = 1.0/pow(loess[i + j*(nrow-1)], 2.0);
			y = -lp2*(e0matrix[i + j*nrow]-m2)/Triangle_c[3 + j*4];
			oneplusexpy = 1+exp(y);
			B = -1.0/oneplusexpy;
			x = -lp1*(e0matrix[i + j*nrow]-m1)/Triangle_c[1  + j*4];
			oneplusexpx = 1+exp(x);
			tmp = k_c[j]*(exp(y)-exp(x))/(oneplusexpy*oneplusexpx);
			A = e0matrix[i + 1 + j*nrow] - e0matrix[i + j*nrow] - tmp;
			a1 += h * pow(B, 2.0);
			b += h * A * B;
			/*c += h * pow(A,2.0);*/
			/*if(j==5)
				Rprintf("\ni=%i: Triangles= %f %f %f %f, e0=%f %f, loess=%f, k=%f", i, 
							Triangle_c[0 + j*4], Triangle_c[1 + j*4], Triangle_c[2 + j*4], Triangle_c[3 + j*4],
							e0matrix[i + j*nrow], e0matrix[i + 1 + j*nrow], 
							loess[i + j*(nrow-1)], k_c[j]);*/
		}
	}
	b = 2*b;
	d[0] = -b/(2*a1);
	d[1] = a1;
}

double rnormtrunc(double mu, double sigma, double low, double high){
	double temp;
	int maxit, i;
	
	temp = -999;
  	maxit = 10;
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


void do_z_gibbs(int *ncountries, int *ntime, double *e0matrix, double *loess, 
						double *Triangle_c, double *k_c, double *p1, double *p2,
						double *omega, double *low, double *high, double *z) {
	double d[2], rn;
	get_z_gibbs_pars(ncountries, ntime, e0matrix, loess, Triangle_c, k_c, p1, p2, d);
	/*Rprintf("\nd = %f, a = %f; sd=%f", d[0], d[1], *omega/pow(2*d[1], 1/2.0));*/
	rn = rnormtrunc(d[0], *omega/pow(2*d[1], 1/2.0), *low, *high);
	/*Rprintf("\nz = %f\n", rn);*/
	z[0] = rn;
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
			double *le, int *idx, int *lidx, double *dct, double *loess_sd, double *logdens) {
	double dl[*lidx], dens[*lidx], dnt[1];
	double s;
	int i;
	dlpars[*par_idx-1] = *x;
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

