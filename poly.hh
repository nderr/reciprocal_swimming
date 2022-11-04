///////////////////////////////////////////////////////////////////////////////
// poly.hh                                                                   //
//                                                                           //
// The poly class represents a set of basis polynomials for a finite element //
//                                                                           //
// Nick Derr                                                                 //
// cleanup 01/29/21                                                          //
///////////////////////////////////////////////////////////////////////////////

#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_poly.h>
#include <cstring>
#include <cstdio>
#include <cmath>

/**
 * Class representing a family of Jacobi polynomials P_n^{a,b} for
 * two provided values a and b, with n in the range [0,max_p]
 *
 * Note that for a=b=-1/2 these are the Chebyshev polynomials and
 * for a=b=0 these are the Legendre polynomials
 */
class jacobi {

	// helper function for applying recurrence relation for
	// the caclulation of polynomial coefficients
	static void jacobi_recur(int n,int a,int b,double *cm2,double *cm1,double *cc);

	public:

	// helper function for calculating the Jacobi coefficients for
	// either a single polynomial or a table of orders [0,p]
	static void jacobi_coeff(int p,double a,double b,double *c,bool keep);

	/** maxmimum polynomial degree */
	int max_p;
	/** Jacobi parameters */
	double alpha,beta;
	/** pointer to coefficient storage */
	double *cf;

	//constructor / destructor
	jacobi(int p,double a,double b);
	~jacobi() {delete[] cf;}

	/** returns a pointer to the coefficients of the pth order polynomial */
	inline double *coeffs(int p) {return cf+p*(max_p+1);}
};

/**
 * Abstract class representing a set of 2D polynomial basis functions
 */
class poly_2d {

	public:
	/** maxmimum polynomial degree */
	int p;
	/** number of functions */
	int n;
	/** pointer to coefficients */
	double *c;
	/** extra space for storage of intermediate coefficient values */
	double *cx,*cdx,*cddx;

	// constructor / destructor
	poly_2d(int p_):p(p_) {};
	virtual ~poly_2d() {};

	// prints out tree of polynomail coefficients
	virtual void print_coeffs(int i)=0;

	// function evaluation
	double eval(int i,double x,double y);
	virtual void eval(int i,double x,double y,double &f,double &fx,double &fy)=0;
	virtual void eval(int i,double x,double y,double &f,double &fx,double &fy,
		double &fxx,double &fxy,double &fyy)=0;

	// polynomial coefficient setup and storage
	virtual void setup_coeffs()=0;
	virtual double* coeffs(int i)=0;
};

/**
 * Class representing a set of Koornwinder polynomials comprising a set
 * of orthonormal polynomials on the standard triangle [0,0],[1,0],[0,1].
 */
class koornwinder : public poly_2d {

	/**
	 * helper function for calculating the non-separable portion of
	 * the (m,k-m) Koornwinder polynomial, i.e.
	 *
	 * P_m^{0,0}(r) * ((1-s)/2)^m
	 *
	 * where r = 2(1+xi)/(1-eta) - 1 and s = eta.
	 *
	 * The coefficients form a flipped upper diagonal matrix,
	 * stored with no zero padding in the array pointed at by cm
	 */
	static void koornwinder_part1(int k,int m,double *cm,double *c);

	public:

	// helper function for calculating the coefficients
	// of the (k,m) Koornwinder polynomial
	static void koornwinder_coeff(int i,jacobi &j00,double *c);
	void setup_coeffs();

	// constructor/destructor
	koornwinder(int p_);
	~koornwinder();

	double *coeffs(int i);
	double *coeffs(int k,int m);

	// prints out tree of polynomail coefficients
	void print_coeffs(int i);


	// functions for polynomial evaluation at provided coordinates
//	double eval(int i,double x,double y);
	void eval(int i,double x,double y,double &f,double &fx,double &fy);
	void eval(int i,double x,double y,double &f,double &fx,double &fy,
		double &fxx,double &fxy,double &fyy);

	static inline void i2km(int i,int &k,int &m) {
		k=static_cast<int>(0.5*(sqrt(1+8*i)-1));
		m=i-k*(k+1)/2;
	}
	static inline int km2i(int k,int m) {return k*(k+1)/2 + m;}
};

/**
 * Class representing a set of Koornwinder polynomials comprising a set
 * of orthonormal polynomials on the standard triangle [0,0],[1,0],[0,1].
 */
class legendre : public poly_2d {

	public:

	// helper function for calculating the coefficients of the (k,m)
	// Koornwinder polynomial
	static void legendre_coeff(int i,jacobi &j00,double *c);
	void setup_coeffs();

	// constructor/destructor
	legendre(int p_);
	~legendre();

	double *coeffs(int i);
	double *coeffs(int k,int m);

	// prints out tree of polynomail coefficients
	void print_coeffs(int i);

	// functions for polynomial evaluation at provided coordinates
//	double eval(int i,double x,double y);
	void eval(int i,double x,double y,double &f,double &fx,double &fy);
	void eval(int i,double x,double y,double &f,double &fx,double &fy,
		double &fxx,double &fxy,double &fyy);

};
