#ifndef MASTER_HH
#define MASTER_HH

#include <cstdio>
#include "poly.hh"
#include "quad.hh"
#include "lapack.hh"
#include "fem.hh"

// forward declarations
struct src_base;
template <typename T> struct tsrc_base;

/**
 * An abstract class representing a 2D master element
 * in a (xi,eta) coordinate system
 */
class master {

	public:
	/** polynomial order */
	int p;
	/** number of nodes */
	int n;
	/** number of gauss quadrature nodes */
	int ng;
	/** pointer to gauss quadrature functions */
	gauss_quad_2d *g;
	poly_2d *polys;

	// matrix storage

	/**
	 * n (rows-i) coefficients of n (cols-j) functions
	 *
	 * ith basis polynomial, jth Lagrange polynomial
	 */
	double *lagr;
	/**
	 * element shape functions
	 *
	 * g (rows-i) functions of n (cols-j) coefficients,
	 * gth Gauss node, jth basis polynomial
	 */
	double *sh;
	/** element shape function xi-derivatives */
	double *sh_xi;
	/** element shape function eta-derivatives */
	double *sh_et;
	/**
	 * element shape functions
	 *
	 * g (rows-i) functions of n (cols-j) coefficients,
	 * gth Gauss node, jth Lagrange polynomial
	 */
	double *sh_L;
	/** element shape function xi-derivatives */
	double *sh_L_xi;
	/** element shape function eta-derivatives */
	double *sh_L_et;

	/**
	 * helper function for printing matrices
	 *
	 * @param[in] m pointer to matrix
	 * @param[in] nr rows in matrix
	 * @param[in] nc cols in matrix
	 * @param[in] out filestream to print to
	 */
	static void print_mat(double *m,int nr,int nc,FILE *out=stdout) {
		for(int j=0;j<nr;j++) {
			for (int i=0;i<nc;i++) {
				fprintf(out," %10.4g ",m[i+j*nc]);
			}
			fprintf(out,"\n");
		}
	}

	// constructors / destructor

	/** constructs a new master element of polynomial order p */
	master(int p_) : p(p_),ng((4*p+3)/2) {}
	/** constructs a new master element of polynomial order p with
	 *  gauss quadrature points corresponding to polynomial order gp */
	master(int p_,int gp) : p(p_),ng((4*gp+3)/2) {}
	virtual ~master();

	/** returns the coordinate of the ith node */
	virtual void point(int i,double &xi,double &eta)=0;
	/** returns the index of the corners */
	virtual void corners(int *crn)=0;
	/** returns the points on the specified side */
	virtual void side(int si,int *sid)=0;
	/** returns the number of sides */
	virtual int n_sides()=0;
	/** returns a pointer to a set of polynomials */
	virtual poly_2d* get_new_poly()=0;
	/** returns a pointer to a gauss quadrature struct */
	virtual gauss_quad_2d* get_new_quad(int ng)=0;

	/** driver method */
	void verify();

	// setup and cleanup functions
	void setup_quad() {g=get_new_quad(ng);}
	void setup_shape_functions();
	void cleanup_shape_functions() {delete[] lagr;delete polys;}
	void cleanup_quad() {delete g;}

	template<typename T>
	void src_term(const trafo_base &ct,const tsrc_base<T> &src,T *s) {

		T zero(0);

		// zero output vector
		for(int i=0;i<n;i++) s[i]=zero;

		// preallocate variables
		double wght,det;
		double xi,et,x,y;
		double sh_pt;
		T src_pt;

		// loop through gauss quad points
		for (int i=0;i<g->N;i++) {

			// grab weighting, determinant, coords and source value
			g->point(i,xi,et);
			wght = g->weight(i);

			// grab Jacobean and coordinate info
			det = ct.Jdet(xi,et);
			ct.xi2x(xi,et,x,y);

			// grab source function strength
			src_pt = -src(x,y);

			T pre=wght*det*src_pt;

			// loop through vector elements
			for(int j=0;j<n;j++) {

				// get contribution of jth Lagrange function at
				// ith quadrature point
				sh_pt = sh_L[j+i*n];

				// add ith quadrature point's contribution
				// to jth vector element
				s[j] += sh_pt*pre;
			}
		}
	}

	/** returns a source term vector */
	void source_term(const trafo_base &ct,const src_base &src,double *s);
	/** returns the element's mass matrix */
	void mass_matrix(const trafo_base &ct,double *m,bool add=false);
	/** returns the element's mass matrix */
	void mass_matrix_cyl_y(const trafo_base &ct,double *m,bool add=false);
	/** returns the element's stiffness matrix */
	void stiffness_matrix(const trafo_base &ct,double *m,bool add=false);
	/** calculates the off-diagonal Bx matrix for a Stokes-type problem */
	static void stokes_bmats(master &mu,master &mp,const trafo_base &ct,
			double *bx,double *by,bool cyl,bool add=false);
	static void check_stokes_up(master &mu,master &mp);

	void L_eval(int j,double xi,double et,double &L,double &L_xi,double &L_et) {

		L = L_xi = L_et = 0;
		double p,p_xi,p_et;

		// fill values
		for(int k=0;k<n;k++) {

			// get poly val
			polys->eval(k,2*xi-1,2*et-1,p,p_xi,p_et);

			// contrib of kth poly to jth lagr
			double M = lagr[j+k*n];

			L    += M*p;
			L_xi += M*p_xi*2; // multiplier because polys \in [-1,1] but using [0,1]
			L_et += M*p_et*2;
		}
	}

	void L_eval(int j,double xi,double et,double &L,double &L_xi,double &L_et,
			double &L_xi_xi,double &L_xi_et,double &L_et_et) {

		L = L_xi = L_et = L_xi_xi = L_xi_et = L_et_et = 0;
		double p,p_xi,p_et,p_xi_xi,p_xi_et,p_et_et;

		// fill values
		for(int k=0;k<n;k++) {

			// get poly val
			polys->eval(k,2*xi-1,2*et-1,p,p_xi,p_et,p_xi_xi,p_xi_et,p_et_et);

			// contrib of kth poly to jth lagr
			double M = lagr[j+k*n];

			L    += M*p;
			L_xi += M*p_xi*2; // multiplier because polys \in [-1,1] but using [0,1]
			L_et += M*p_et*2;

			L_xi_xi += M*p_xi_xi*4;
			L_xi_et += M*p_xi_et*4;
			L_et_et += M*p_et_et*4;
		}
	}

	template <typename T>
	void eval(trafo_base &ct,double xi,double et,T *dat,T &v,T &v_x,T &v_y) {

		// out vals + Lagrange poly vals
		v = v_x = v_y = 0;
		double L,L_xi,L_et;

		// inv Jacobeans
		double xi_x,et_x,xi_y,et_y;
		ct.Jinv(xi,et,xi_x,et_x,xi_y,et_y);

		for (int l=0;l<n;l++) {

			// grab data
			T d = dat[l];

			// grab Lagrange vals
			L_eval(l,xi,et,L,L_xi,L_et);

			// linear combination, weighted
			// by Lagrange vals
			v   +=                       L*d;
			v_x += (L_xi*xi_x + L_et*et_x)*d;
			v_y += (L_et*et_y + L_xi*xi_y)*d;
		}
	}

	template <typename T>
	void eval(trafo_base &ct,double xi,double et,T *dat,T &v,T &v_x,T &v_y,
			T &v_xx,T &v_xy,T &v_yy,bool pr=false) {

//		if (pr) for(int i=0;i<n;i++) printf("v[%d]=%g + i*%g\n",i,std::real(dat[i]),std::imag(dat[i]));

		// out vals + Lagrange poly vals
		v = v_x = v_y = v_xx = v_xy = v_yy = 0;
		double L,L_xi,L_et,L_xi_xi,L_xi_et,L_et_et;

		// inv Jacobeans
		double xi_x,et_x,xi_y,et_y;
		ct.Jinv(xi,et,xi_x,et_x,xi_y,et_y);

		double xi_xx,xi_xy,xi_yy,et_xx,et_xy,et_yy;
		ct.Jinv2(xi,et, xi_xx, xi_xy, xi_yy, et_xx, et_xy, et_yy);

		for (int l=0;l<n;l++) {

			// grab data
			T d = dat[l];

			// grab Lagrange vals
			L_eval(l,xi,et,L,L_xi,L_et,L_xi_xi,L_xi_et,L_et_et);

			// linear combination, weighted
			// by Lagrange vals
			v   +=                       L*d;
			v_x += (L_xi*xi_x + L_et*et_x)*d;
			v_y += (L_et*et_y + L_xi*xi_y)*d;

//
			v_xx += (L_xi * xi_xx + L_et * et_xx + L_xi_xi * xi_x * xi_x
					+ L_et_et * et_x * et_x + 2 * L_xi_et * xi_x * et_x) * d;

			v_xy += (L_xi * xi_xy + L_et * et_xy + L_xi_xi * xi_x * xi_y
					+ L_et_et * et_x * et_y + L_xi_et * (xi_x * et_y + xi_y * et_x)) * d;

			v_yy += (L_xi * xi_yy + L_et * et_yy + L_xi_xi * xi_y * xi_y
					+ L_et_et * et_y * et_y + 2 * L_xi_et * xi_y * et_y) * d;

//.			if (pr) printf("cxx[%d]=%g\n",l,
//(L_xi * xi_xx + L_et * et_xx + L_xi_xi * xi_x * xi_x
//					+ L_et_et * et_x * et_x + 2 * L_xi_et * xi_x * et_x)
		}
	}

	template <typename T>
	void eval(double xi,double et,T *dat,T &v,T &v_xi,T &v_et) {

		v = v_xi = v_et = 0;
		double L,L_xi,L_et;

		for (int l=0;l<n;l++) {

			// grab data
			T d = dat[l];

			// grab Lagrange vals
			L_eval(l,xi,et,L,L_xi,L_et);

			// linear combination, weighted
			// by Lagrange vals
			v    +=    L*d;
			v_xi += L_xi*d;
			v_et += L_et*d;
		}
	}

	/** evalutes solution at a given point */
	template <typename T>
	T eval(double xi,double et,T *dat) {
		
		// get polynomial values at this point
		// TODO just put canned function in polys calss
		// for direct access to Lagrange...would be easier to read
		double *poly_vals = new double[n];
		for(int i=0;i<n;i++) poly_vals[i] = polys->eval(i,2*xi-1,2*et-1);

		// now convert to Lagrange basis
		double *lagr_vals = new double[n];
		for(int i=0;i<n;i++) lagr_vals[i]=0;

		// contribution of kth poly to jth lagr
		for(int k=0;k<n;k++)for(int j=0;j<n;j++)
			lagr_vals[j] += lagr[j+k*n]*poly_vals[k];

		// weighted average of provided data points
		T out(0);
		for(int i=0;i<n;i++) out += lagr_vals[i] * dat[i];

		delete[] lagr_vals;
		delete[] poly_vals;

		return out;
	}

	/** prints out element boundaries */
	void gp_lines(const trafo_base &ct,double *dat,FILE *out=stdout);
	/** prints out scalar function value over domain */
	void gp_fill(const trafo_base &ct,double *dat,FILE *out=stdout);
};

/**
 * Class representing a 2D quad master element
 */
class master_quad : public master {
	public:
	/** number of sides */
	static const int ns=4;

	// constructors
	master_quad(int p,int ng):master(p,ng){
		setup_quad();
		setup_shape_functions();
	}
	master_quad(int p):master(p){
		setup_quad();
		setup_shape_functions();
	}


	/** returns the coordinate of the ith node */
	void point(int i,double &xi,double &eta);
	/** returns the index of the corners */
	void corners(int *crn);
	/** returns the points on the specified side */
	void side(int si,int *sid);
	/** returns the number of sides */
	int n_sides() {return ns;}
	/** returns a pointer to a set of polynomials */
	poly_2d* get_new_poly() {return new legendre(p);}
	/** returns a pointer to a gauss quadrature struct */
	gauss_quad_2d* get_new_quad(int gp) {return new gauss_quad_2d_square_tens(gp);}

	/** 
	 * calculates interpolation matrix from one element to 4
	 * @param[out] i_mat interpolation matrix. (2p+1)^2 rows,
	 *   representing points in the fine element and
	 *   (p+1)^2 cols, representing points in coarse element
	 */
	void alloc_interp(double *(&i_mat)) {

		// polynomial values
		double *p_vals = new double[n];

		// interpolation values
		i_mat = new double[(p+1)*(p+1) * (2*p+1)*(2*p+1)];

		// for each fine point 
		for (int fr=0,row=0;fr<=2*p;fr++) for (int fc=0;fc<=2*p;fc++,row++) {

			// fine local points
			double xi_f=(1.*fc)/(2*p), eta_f=(1.*fr)/(2*p);

			// get each poly basis val at this point
			for(int i=0;i<n;i++)
				p_vals[i] = polys->eval(i,2*xi_f-1,2*eta_f-1);

			// convert to shape function representation
			double *l_vals = i_mat + row*n;  // points to one row of conv mat
			for(int i=0;i<n;i++) l_vals[i]=0;

			// (jth Lagrange poly, kth poly basis function)
			for (int j=0;j<n;j++) for (int k=0;k<n;k++)
				l_vals[j] += lagr[j+k*n]*p_vals[k];

			// regularize
			for(int i=0;i<n;i++) if (fabs(l_vals[i])<1e-5) l_vals[i]=0;
		}

		// now, i_mat[j + k*n] is the jth coarse point's
		// contribution to the kth fine point

		delete[] p_vals;
	}

	void free_interp(double *(&i_mat)) {delete[] i_mat; i_mat=NULL;}

	/**
	 * Returns the matrices \int \psi_i\psi_k\partial_x(\psi_j) dA and
	 *   \int \psi_i\psi_k\partial_y(\psi_j) dA.
	 *
	 * Note this is an n x n x n matrix
	 *
	 * @param[in] ct coordinate transform for the element
	 * @param[out] mx the matrix corresponding to x-derivs of the shape functions
	 * @param[out] my the matrix corresponding to y-derivs of the shape functions
	 */
	void reciprocal_thm_mat(trafo_base &ct,double *mx,double *my,bool add=false) {

		// if not adding in, zero out values
		if (!add) for(int i=0;i<n*n*n;i++) {mx[i]=my[i]=0;}

		// preallocate variables
		double wght,det;
		double xi,et;
		double xi_x,xi_y,et_x,et_y;
		double sh_i,sh_k,sh_x_j,sh_y_j,sh_xi,sh_et;

		// march through gauss nodes
		for (int gi=0;gi<g->N;gi++) {

			// grab coord, weight and determinant
			g->point(gi,xi,et);
			wght = g->weight(gi);

			// get determinant and inverse Jacobean
			det = ct.Jdet(xi,et);
			ct.Jinv(xi,et,xi_x,et_x,xi_y,et_y);

			// march through i nodes
			for(int i=0;i<n;i++) {

				sh_i = sh_L[i + gi*n];

				// march through k nodes
				for(int k=0;k<n;k++) {

					sh_k = sh_L[k + gi*n];

					// march through j nodes
					for(int j=0;j<n;j++) {

						// grab local derivs
						sh_xi = sh_L_xi[j+gi*n];
						sh_et = sh_L_et[j+gi*n];

						// convert to global derivs
						sh_x_j = xi_x * sh_xi + et_x * sh_et;
						sh_y_j = xi_y * sh_xi + et_y * sh_et;

						// add into matrix
						mx[i + j*n + k*n*n] += sh_i * sh_k * sh_x_j * det * wght;
						my[i + j*n + k*n*n] += sh_i * sh_k * sh_y_j * det * wght;
					}
				}
			}
		}
	}
};

/**
 * Class representing a 2D quad master element
 */
class master_tri : public master {
	public:
	/** number of sides */
	static const int ns=3;
	
	// constructor
	master_tri(int p):master(p){
		setup_quad();
		setup_shape_functions();
	}

	/** returns the coordinate of the ith node */
	void point(int i,double &xi,double &eta);
	/** returns the index of the corners */
	void corners(int *crn);
	/** returns the points on the specified side */
	void side(int si,int *sid);
	/** returns the number of sides */
	int n_sides() {return ns;}
	/** returns a pointer to a set of polynomials */
	poly_2d* get_new_poly() {return new koornwinder(p);}
	/** returns a pointer to a gauss quadrature struct */
	gauss_quad_2d* get_new_quad(int gp) {return new gauss_quad_2d_square_tens(gp);}

};

/**
 * Abstract class for a functor representing a position-dependent
 * source function
 */
struct src_base {
	virtual double operator()(double x,double y)const=0;
};

template <typename T>
struct tsrc_base {
	virtual T operator()(double x,double y)const=0;
};

#endif
