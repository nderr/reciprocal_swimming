///////////////////////////////////////////////////////////////////////////////
// poly.cc                                                                   //
//                                                                           //
// Implementations for the  poly class represents a set of basis polynomials //
// for a finite element                                                      //
//                                                                           //
// Nick Derr                                                                 //
// cleanup 01/29/21                                                          //
///////////////////////////////////////////////////////////////////////////////

#include "poly.hh"

jacobi::jacobi(int p,double a,double b) : max_p(p),alpha(a),beta(b),
		cf(new double[(max_p+1)*(max_p+1)]) {
		jacobi_coeff(p,a,b,cf,true);
}

double* legendre::coeffs(int i) {
	int j,k;
	j=i/(p+1);
	k=i%(p+1);
	return coeffs(j,k);
}

double* legendre::coeffs(int j,int k) {
	double *cc=c;
	for(int jj=0;jj<=j;jj++) {
		for(int kk=0;kk<=(j==jj?k-1:p);kk++) {
			cc+=(jj+1)*(kk+1);
		}
	}
	return cc;
}

double* koornwinder::coeffs(int i) {
	int k,m;
	i2km(i,k,m);
	return coeffs(k,m);
}

double* koornwinder::coeffs(int k,int m) {
	// dist to start of k region
	int to_k = k*(1+k)*(1+k)*(2+k)/12;
	// dist to m within this k region
	int to_m = (4+3*k-m)*m*(1+m)/6;
	return c+to_k+to_m;
}

legendre::legendre(int p) : poly_2d(p) {setup_coeffs();}

koornwinder::koornwinder(int p) : poly_2d(p) {setup_coeffs();}

void legendre::setup_coeffs() {

	n = (p+1)*(p+1);
	c = new double[(1+p)*(1+p)*(2+p)*(2+p)/4];
	cx = new double[p+1];
	cdx = new double[p+1];
	cddx = new double[p+1];

	jacobi j(p,0,0);
	for(int i=0;i<n;i++) legendre_coeff(i,j,coeffs(i));
}

void koornwinder::setup_coeffs() {

	n = (p+1)*(p+2)/2;
	c = new double[(1+p)*(2+p)*(2+p)*(3+p)/12];
	cx = new double[p+1];
	cdx = new double[p+1];
	cddx = new double[p+1];

	jacobi j(p,0,0);
	for(int i=0;i<n;i++)koornwinder_coeff(i,j,coeffs(i));
}

koornwinder::~koornwinder() {
	delete[] cddx;
	delete[] cdx;
	delete[] cx;
	delete[] c;
}

legendre::~legendre() {
	delete[] cddx;
	delete[] cdx;
	delete[] cx;
	delete[] c;
}

/** helper function for printing out polynomial coefficients */
void legendre::print_coeffs(int i) {
	double *cc = coeffs(i);
	int j,k;
	j=i/(p+1);
	k=i%(p+1);

	printf("%d:\n",i);
	for (int jj=0;jj<=j;jj++) {
		for (int kk=0;kk<=k;kk++,cc++) {
			printf(" %g ",*cc);
		}
		printf("\n");
	}
	printf("\n");
}

/** helper function for printing out polynomial coefficients */
void koornwinder::print_coeffs(int i) {
	double *cc = coeffs(i);
	int k,m;
	i2km(i,k,m);
	printf("%d:\n",i);
	for (int l=0;l<=m;l++) {
		for (int p=0;p<=k-l;p++,cc++) {
			printf(" %g ",*cc);
		}
		printf("\n");
	}
	printf("\n");
}

double poly_2d::eval(int i,double x,double y) {
	double f,trash1,trash2;
	eval(i,x,y,f,trash1,trash2);
	return f;
}


/**
 * evalutes the i-th polynomial at the point (x,y), calculating the
 * function value f and the gradient components (fx,fy).
 */
void koornwinder::eval(int i,double x,double y,double &f,double &fx,double &fy) {

	double res[2];

	int k,m;
	i2km(i,k,m);
	double *cc = coeffs(k,m);

	for (int j=0;j<=m;cc+=(k-j+1),j++) {
		gsl_poly_eval_derivs(cc,k-j+1,y,res,2);
		cx[j]=res[0];
		cdx[j]=res[1];
	}
	
	gsl_poly_eval_derivs(cx,m+1,x,res,2);
	f = res[0];
	fx = res[1];
	fy = gsl_poly_eval(cdx,m+1,x);
}

/**
 * evalutes the i-th polynomial at the point (x,y), calculating the
 * function value f and the gradient components (fx,fy).
 */
void koornwinder::eval(int i,double x,double y,double &f,double &fx,double &fy,
	double &fxx,double &fxy,double &fyy) {

	double res[3];

	int k,m;
	i2km(i,k,m);
	double *cc = coeffs(k,m);

	for (int j=0;j<=m;cc+=(k-j+1),j++) {
		gsl_poly_eval_derivs(cc,k-j+1,y,res,3);
		cx[j]=res[0];
		cdx[j]=res[1];
		cddx[j]=res[2];
	}
	
	gsl_poly_eval_derivs(cx,m+1,x,res,3);
	f = res[0];
	fx = res[1];
	fxx = res[2];

	gsl_poly_eval_derivs(cdx,m+1,x,res,2);
	fy = res[0];
	fxy = res[1];

	fyy = gsl_poly_eval(cddx,m+1,x);
}

void legendre::eval(int i,double x,double y,double &f,double &fx,double &fy) {

	double res[2];

	int j=i/(p+1);
	int k=i%(p+1);

	double *cc = coeffs(i);
	for (int jj=0;jj<=j;cc+=(k+1),jj++) {
		gsl_poly_eval_derivs(cc,k+1,y,res,2);
		cx[jj]=res[0];
		cdx[jj]=res[1];
	}
	
	gsl_poly_eval_derivs(cx,j+1,x,res,2);
	f = res[0];
	fx = res[1];
	fy = gsl_poly_eval(cdx,j+1,x);
}

void legendre::eval(int i,double x,double y,double &f,double &fx,double &fy,
	double &fxx,double &fxy,double &fyy) {

	double res[3];

	int j=i/(p+1);
	int k=i%(p+1);

	double *cc = coeffs(i);
	for (int jj=0;jj<=j;cc+=(k+1),jj++) {
		gsl_poly_eval_derivs(cc,k+1,y,res,3);
		cx[jj]=res[0];
		cdx[jj]=res[1];
		cddx[jj] = res[2];
	}
	
	gsl_poly_eval_derivs(cx,j+1,x,res,3);
	f = res[0];
	fx = res[1];
	fxx = res[2];

	gsl_poly_eval_derivs(cdx,j+1,x,res,2);
	fy = res[0];
	fxy = res[1];

	fyy = gsl_poly_eval(cddx,j+1,x);
}

/**  
 *   applies the recursion relation
 *   
 *   an*f_n = (a1*x+a2)*f_{n-1} + a3*f_{n-2}
 *
 *   where f_n is the Jacobi polynomial P_n^{(a,b)}
 *
 *   to obtain the nth degree polynomial coefficients
 */
void jacobi::jacobi_recur(int n,int a,int b,double *cm2,double *cm1,double *cc) {

	// "base" values
	const double b1 = n+a+b-1;
	const double b2 = b1+n-1;
	const double b3 = a*a-b*b;

	// use "base" values to construct factors from recursion relation
	const double an = 2*n*(b1+1)*b2;
	const double a1 = b2*(b2+1)*(b2+2);
	const double a2 = (b2+1)*b3;
	const double a3 = -2*(b1-b)*(b1-a)*(b2+2);
	
	double v1,v2,v3;
	for(int i=0;i<=n;i++){

		// bit from x * f_{n-1}
		v1 = i>0?a1*cm1[i-1]:0;

		// bit from f_{n-1}
		v2 = i==n?0:a2*cm1[i];

		// bit from f_{n-2}
		v3 = i<n?a3*cm2[i]:0;

		// calculate f_n
		cc[i] = (v1+v2+v3)/an; 
	}
}

/** calculates the polynomial coefficients of the nth-degree
 *  Jacobi polynomial P_n^{a,b}. Optionally returns the coefficients
 *  for each polynomial of degree less than or equal to n */
void jacobi::jacobi_coeff(int n,double a,double b,double *c,bool keep=false) {

	// number of terms in polynomial 
	const int nt = n+1;

	// base cases
	const double cm20 = 1;
	const double cm10 = 0.5*(a-b);
	const double cm11 = 1+0.5*(a+b);

	// check base cases
	if (n==0) {

		// if only one term, doesn't matter if keeping or not
		*c = cm20;
		return;

	} else if (n==1) {

		// if linear case, return matr or vec as requested
		if (keep) {
			*c = cm20;
			c[1] = 0;
			c[2] = cm10;
			c[3] = cm11;
		} else {
			*c = cm10;
			c[1] = cm11;
		}
		return;
	}

	// pointers to n-2, n-1 and n coefficients
	// only allocate if not returning all coeffs
	double *cm2 = keep ?  c      : new double[nt],
		   *cm1 = keep ? (c+nt)  : new double[nt],
		   *cc  = keep ? (cm1+nt): c;

	// fill in with initial cases
	for(int i=0;i<nt;i++) {
		cm2[i] = i==0?cm20:0;
		cm1[i] = i==0?cm10:(i==1?cm11:0);
	}

	// calculate each order of poly coeffs in turn
	for (int lev=2;lev<nt;lev++) {

		// recur and fill in remainder with zeros
		jacobi_recur(lev,a,b,cm2,cm1,cc);
		for(int i=lev+1;i<nt;i++) cc[i]=0;

		if (keep) {

			// if keeping, move through array
			cm2=cm1;
			cm1=cc;
			cc+=nt;
		} else {

			// otherwise, swap pointers and copy in
			double *tmp=cm2;
			cm2=cm1;
			cm1=tmp;	
			memcpy(cm1,cc,nt*sizeof(double));
		}
	}

	// clean up if not returning everything
	if (!keep) {
		delete[] cm1;
		delete[] cm2;
	}
}


void koornwinder::koornwinder_part1(int k,int m,double *cm,double *c) {

	double mpre = pow(0.5,m),ppre,jpre,sgn;

	for(int l=0;l<=m;l++) {
		for(int p=0;p<=m-l;p++,c++) {
			*c=0;
			ppre = p%2==0?1:-1;
			for(int j=l;j<=m-p;j++) {
				jpre = (1<<j)*gsl_sf_choose(m-j,p)*gsl_sf_choose(j,l);
				for (int i=j;i<=m;i++) {
					sgn = ((i-j)%2==0)?1:-1;
					*c += mpre*ppre*jpre*cm[i]*sgn*gsl_sf_choose(i,j);
				}
			}
		}
	}
}

void legendre::legendre_coeff(int ip,jacobi &j00,double *c) {

	int j,k;
	j = ip/(j00.max_p+1);
	k = ip%(j00.max_p+1);

	// normalization constant
	double norm = sqrt((2*j+1)*(2*k+1));
	double *cj = j00.coeffs(j);
	double *ck = j00.coeffs(k);

	// each x coord (each row)
	for (int jj=0;jj<=j;jj++) {
		for (int kk=0;kk<=k;kk++,c++) {
			*c = cj[jj]*ck[kk]*norm;
		}
	}
}

void koornwinder::koornwinder_coeff(int ip,jacobi &j00,double *c) {

	int k,m;
	i2km(ip,k,m);

	// normalization constant
	double norm = sqrt(2*(2*m+1)*(k+1));

	// grab first part
	double *c1 = new double[(m+1)*(m+2)/2],*cp1=c1;
	koornwinder_part1(k,m,j00.coeffs(m),c1);

	// now convolve
	double *cn = new double[k-m+1];
	jacobi::jacobi_coeff(k-m,2*m+1,0,cn);

	// each x coord (each row)
	for (int i=0;i<=m;cp1+=(m-i+1),i++) {
		for (int j=0;j<=k-i;j++,c++) {
			*c = 0;
			for (int jc=m-k;jc<=0;jc++) {
				int ii = j+jc;
				*c += (ii>=0&&ii<=m-i)?cp1[ii]*cn[-jc]:0;
			}
			*c *= norm;
		}
	}

	delete[] cn;
	delete[] c1;
}
