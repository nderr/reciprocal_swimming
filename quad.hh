///////////////////////////////////////////////////////////////////////////////
// quad.hh
//
// collection of quadrature rules
///////////////////////////////////////////////////////////////////////////////
#ifndef QUAD_HH
#define QUAD_HH

#include <gsl/gsl_integration.h>
#include <cstdio>
#include <cstring>

/** class representing a 1D Gauss-Legendre fixed-point
 *  quadrature rule on the interval [0,1]. An n-point
 *  rule will exactly integrate polynomials of degree
 *  less than 2*n - 1 */
class gauss_quad_1d {

	public:
	
	/** number of points */
	int n;
	/** abcissae and weights of the n points */
	double *x,*w;

	// constructor/destructor
	gauss_quad_1d(int n_);
	~gauss_quad_1d();
};

/** abstract class representing a 2D Gauss-Legendre quadrature
 *  rule on some domain */
class gauss_quad_2d {
	public:

	/** total number of points in the quadrature rule */
	int N;
	/** abscissae of quadrature points */
	double *x,*y;
	/** weights of quadrature points */
	double *w;

	// constructors/destructor
	gauss_quad_2d(int N_);
	gauss_quad_2d(int N_,int nx,int ny);
	virtual ~gauss_quad_2d();

	/** returns the weight of the i-th quadrature point */
	inline double weight(int i) {return w[i];};
	/** provides the x and y coordinates of the i-th quadrature point */
	virtual void point(int i,double &xx,double &yy) {xx=x[i];yy=y[i];};
	/** sets up abscissae and weights of quadrature rule (pure) */
	virtual void setup_rule() = 0;
};

/** class representing a 2D quadrature rule over the domain [0,1] x [0,1],
 *  constructed as a tensor product of 1D quadrature rules */
class gauss_quad_2d_square_tens : public gauss_quad_2d {

	/** number of quadrature points in one direction */
	int n;

	public:
	/** constructor */
	gauss_quad_2d_square_tens(int n_): gauss_quad_2d(n_*n_,n_,n_),
		n(n_) {setup_rule();}

	/** constructs 2D quadrature rule as a tensor product of 1D rule */
	void setup_rule();

	/** provides the x and y components of the i-th quadrature point */
	void point(int i,double &xx,double &yy) {xx=x[i%n];yy=y[i/n];}
};

/** class representing a 2D quadrature rule over the standard triangle 
 *  [0,0] -> [1,0] -> [0,1], constructed by mapping the tensor-product
 *  rule on the square [0,1]x[0,1] to the triangular domain */
class gauss_quad_2d_tri_tens : public gauss_quad_2d {

	/** number of quadrature points in one direction */
	int n;

	public:
	/** constructor */
	gauss_quad_2d_tri_tens(int n_) : gauss_quad_2d(n_*n_,n_*n_,n_),
		n(n_) {setup_rule();}

	/** constructs the quadrature rule by mapping the abscissae and
	 *  weights from the square [0,1]x[0,1] to the standard triangle */
	void setup_rule();

	/** provides the x and y components of the i-th quadrature point */
	void point(int i,double &xx,double &yy) {xx=x[i];yy=y[i/n];}
};

#endif
