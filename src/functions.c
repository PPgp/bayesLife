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
