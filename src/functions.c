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

void get_z_gibbs_pars(int *it, int *ncountries, int *ntime, double *e0matrix, double *loess, 
						double *Triangle_c, double *k_c, double *p1, double *p2,
						double *omega, double *alpha, double *delta, double *d) {
	int i, j, nrow, ncol, xisinf, yisinf, l, create_file;
	double A, B, y, x, m1, m2, lp1, lp2, b, c, h, oneplusexpy, oneplusexpx, gterm, a1, xl, yl;
	double thisa1, thisb, thisd, al, del, a, om, ado, delinv, zmean;
	double At[11], Bt[11], ht[11], xlt[11], ylt[11];
	ncol = *ncountries;
	char f[20];
	FILE * pFile;
	a1 = 0;
	b = 0;
	al = *alpha;
	del = *delta;
	om = *omega;
	/*ado = a/pow(om, 2.0);
	delinv = 1.0/pow(del, 2.0);*/
	/*c = 0;*/
	lp1 = log(pow(*p1,2.0));
	lp2 = log(pow(*p2,2.0));
	create_file=0;
	/*if(*it==874 || *it==875 || *it==800 || *it==950) create_file=1;*/
	if(*it==30 || *it==99) create_file=1;
	create_file=1;
	if(create_file) {
		sprintf(f, "AB_%i.txt", *it);
		pFile = fopen(f,"w");
	}
	for (j=0; j<ncol; j++) {
		nrow = ntime[j];
		m1 = Triangle_c[j*4] + 0.5*Triangle_c[1 + j*4];
		m2 = Triangle_c[j*4] + Triangle_c[1 + j*4] + Triangle_c[2 + j*4] + 0.5*Triangle_c[3 + j*4];
		if(create_file)Rprintf("\nit= %i, j=%i", *it, j);
		thisa1 = 0;
		thisb = 0;
		for (i=0; i<(nrow-1); i++) {
			h = 1.0/pow(loess[i + j*(nrow-1)], 2.0);
			xl = -lp1*(e0matrix[i + j*nrow]-m1)/Triangle_c[1 + j*4];
			xisinf = 0;
			if(xl < -700) x = 0;
			else{
				if(xl > 700) xisinf = 1;
				else x = exp(xl);
			}
			yl = -lp2*(e0matrix[i + j*nrow]-m2)/Triangle_c[3 + j*4];
			yisinf = 0;
			if(yl < -700) y = 0;
			else{
				if(yl > 700) yisinf = 1;
				else y = exp(yl);
			}
			if (xisinf || yisinf) gterm = 0;
			else gterm = k_c[j]*(y-x)/((1+x)*(1+y));
			A = e0matrix[i + 1 + j*nrow] - e0matrix[i + j*nrow] - gterm;
			if(yisinf) B=0; else B = -1.0/(1+y);
			thisa1 += h * pow(B, 2.0);
			thisb += h * A * B;
			At[i] = A;
			Bt[i] = B;
			ht[i] = h;
			xlt[i] = xl;
			ylt[i] = yl;
			/*c += h * pow(A,2.0);*/
			 if(create_file) {
			 	Rprintf("\ni=%i: Triangles= %f %f %f %f, e0=%f %f, loess=%f, k=%f", i, 
							Triangle_c[0 + j*4], Triangle_c[1 + j*4], Triangle_c[2 + j*4], Triangle_c[3 + j*4],
							e0matrix[i + j*nrow], e0matrix[i + 1 + j*nrow], 
							loess[i + j*(nrow-1)], k_c[j]);
				Rprintf("\n     h=%f, A=%f, B=%f, a=%f, b=%f", h, A, B, h * pow(B, 2.0), 2*h * A * B);
			}
		}
		/*thisd = -2*thisb/(2*thisa1);
		zmean = (ado*thisd + delinv*al)/(ado + delinv);
		if(zmean >= 0 && zmean <= 1.15) {
			if(create_file) Rprintf("\nACCEPTED");*/
			a1 += thisa1;
			b += thisb;	
			if(create_file) {
				for (l=0; l<(nrow-1); l++) {
					fprintf(pFile, "%i %i %f %f %f %f %f %f %f %f %f %f %f %f\n", l, j, At[l], Bt[l], ht[l], Triangle_c[0 + j*4], Triangle_c[1 + j*4], Triangle_c[2 + j*4], Triangle_c[3 + j*4], k_c[j], e0matrix[l + j*nrow], loess[l + j*(nrow-1)], xlt[l], ylt[l]);
					/*fprintf(pFile, "%i %i %f %f %f\n", i, *it, A, B, h);*/
				}
			}
		/*} else if(create_file) Rprintf("\nNOT ACCEPTED");*/
	}
	if(create_file) fclose(pFile);
	b = 2*b;
	d[0] = -b/(2*a1);
	d[1] = a1;
	d[2] = b;
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


void do_z_gibbs(int *it, int *ncountries, int *ntime, double *e0matrix, double *loess, 
						double *Triangle_c, double *k_c, double *p1, double *p2,
						double *omega, double *alpha, double *delta, 
						double *low, double *high, double *zmean, double *sd_d, double *z) {
	double d[3], rn, al, del, a, om, ado, delinv;
	get_z_gibbs_pars(it, ncountries, ntime, e0matrix, loess, Triangle_c, k_c, p1, p2, 
					omega, alpha, delta, d);

	/*z[0] = rn;*/
	a = d[1];
	al = *alpha;
	del = *delta;
	om = *omega;
	ado = a/pow(om, 2.0);
	delinv = 1.0/pow(del, 2.0);
	zmean[0] = (ado*d[0] + delinv*al)/(ado + delinv);
	sd_d[0] = pow(1.0/(ado + delinv), 1/2.0);
	z[0] = rnormtrunc(zmean[0], sd_d[0], *low, *high);
	/*if(*it==30 || *it==99) {*/
		Rprintf("\nd = %f, a = %f; b = %f; omega=%f; zmean=%fd sd=%f; z=%f", d[0], d[1], d[2], *omega, 
				zmean[0], sd_d[0], z[0]);/*}*/
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
