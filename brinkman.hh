#ifndef BRINKMAN_HH
#define BRINKMAN_HH

#include <cassert>
#include <omp.h>
#include <complex>
#include <vector>
#include <numeric>
#include <limits>

#include "quad.hh"
#include "mpi.h"
#include "fem.hh"
#include "master.hh"
#include "lapack.hh"
#include "petscksp.h"
#include "petscdm.h"
#include "petscdmda.h"
#include "petscdmcomposite.h"

/** helper function - am I on rank 0? */
bool gru();

/** helper function for bailing */
void bail(int code,const char *fmt,...);

/** helper function for printing */
void mesg(const char *fmt,...);

/** complex doubles */
typedef std::complex<double> dcomplex;

/** problem types implemented */
enum class ptype {two_spheres,one_sphere,driven_cyl};

/** vector of complex doubles */
struct vcomplex {

	/** x,y components */
	dcomplex x,y;

	// constructors
	vcomplex():x(0),y(0) {}
	vcomplex(dcomplex v):x(v),y(v) {

		if (std::norm(v) != 0) {
			fputs("bad!\n",stderr);
			bail(-1,"don't use this");
		}
	}
	vcomplex(dcomplex x_,dcomplex y_):x(x_),y(y_) {}

	/** negation operator */
	vcomplex operator-() const {return vcomplex(-x,-y);}
	/** addition operator */
	vcomplex operator+(const vcomplex &v) const {return vcomplex(x+v.x,y+v.y);}
	/** subtraction operator */
	vcomplex operator-(const vcomplex &v) const {return vcomplex(x+(-v.x),y+(-v.y));}
	/** multiplication operator */
	vcomplex operator*(const double &s) const {return vcomplex(s*x,s*y);}
	/** multiplication operator */
	vcomplex operator*(const dcomplex &s) const {return vcomplex(s*x,s*y);}
	/** addition update */
	vcomplex& operator+=(const vcomplex &v) {x+=v.x; y+=v.y; return *this;}
	/** subtraction update */
	vcomplex& operator-=(const vcomplex &v) {return (*this)+=(-v);}
	/** div update */
	vcomplex& operator/=(const dcomplex &s) {x/=s; y/=s; return *this;}
	/** real part */
	vcomplex real() const {return vcomplex(std::real(x),std::real(y));}
	/** imag part */
	vcomplex imag() const {return vcomplex(std::imag(x),std::imag(y));}

	static vcomplex zero;
};


/** commutative multiplication */
static vcomplex operator*(const double &s,const vcomplex &v) {return v*s;}
/** commutative multiplication */
static vcomplex operator*(const dcomplex &s,const vcomplex &v) {return v*s;}

/** imaginary number */
static const dcomplex I = sqrt(dcomplex(-1));

/** homogeneous complex source function */
struct homog:public tsrc_base<vcomplex> {
	vcomplex operator()(double x,double y) const {return vcomplex::zero;}
};

// forward declarations
struct brinkman;
struct petsc_solve;

/** source function for Reynolds stress given a Brinkman solution */
struct reynolds:public tsrc_base<vcomplex> {

	/** Brinkman solution */
	brinkman *br;
	// constructor
	reynolds(brinkman *br_):br(br_) {}
	/** source at location (x,y) */
	vcomplex operator()(double x,double y) const;
};

/** source function for method of manufactured solutions */
/*
struct mmr: public tsrc_base<vcomplex> {

	vcomplex operator()(double x,double y) const {

		vcomplex v(

		return 
	}
};
*/

/**
 * represents the necessary pieces of a Petsc solve, including
 * DMs representing the FE grid, velocity and pressure node grids,
 * and a DMComposite representing total nodal grid. Also has the
 * problem matrix, source and solution vectors, and solver
 */
struct petsc_solve {

	/** wrapper for brinkman problem description */
	struct user_context {
		brinkman *prob;
	} ctx;
	/** linear solver */
	KSP ksp;

	/** solution method type */
	enum class method {lu,schur,mg};

	// constructor
	petsc_solve(brinkman *prob,method type=method::lu);
	~petsc_solve() { KSPDestroy(&ksp);}

	// helpers
	
	// allocate system DM
	static void alloc_sys_dmda(brinkman*,DM&);
	// store user context for operator construction
	void assign(brinkman *prob) {ctx.prob=prob;}
	// apply Petsc command-line args
	void add_opt(const char* s1,const char *s2) {
		PetscOptionsSetValue(NULL,s1,s2);}
	void add_opt(const char* s1) {
		PetscOptionsSetValue(NULL,s1,NULL);}
	
	// set up methods with Petsc options
	void set_lu(const char* lu_type="mumps");
	void set_schur(const char* lu_type="mumps");
	void set_mg(const char* lu_type="mumps");

	// solve the actual brinkman problem
	void solve();

	// calculate the pressure node mass matrix
	static PetscErrorCode calculate_p_mass(KSP ksp,Mat &P,void *ctx);
	// calculate the pressure node stiffness matrix
	static PetscErrorCode calculate_p_stiff(KSP ksp,Mat &P,void *ctx);
	// calculate the problem matrix (composite or vel)
	static PetscErrorCode calculate_A(KSP ksp,Mat A,Mat P,void *ctx);

	// calculate the velocity part of the problem matrix
	static void calculate_A_vel(DM da_vel,Mat P,brinkman*);
	// calculate the pressure part of the problem matrix
//	static void calculate_A_pre(DM da_pre,Mat P,brinkman*);
	// calculate the total problem matrix
	static void calculate_A_sys(DM da_sys,Mat P,brinkman*);

	// caclulate the RHS vector
	static PetscErrorCode calculate_B(KSP ksp,Vec B,void *ctx);
//	static PetscErrorCode calculate_AU(KSP ksp,Mat A,Mat P,void *ctx);
//	static PetscErrorCode calculate_BU(KSP ksp,Vec B,void *ctx);

	void check_sol() {

		KSPConvergedReason reason;
		KSPGetConvergedReason(ksp,&reason);

		if (reason>0) return;

		
		
		KSPDestroy(&ksp);
		bail(-1,"didn't converge for reason %d\n",reason);
	}

	static int sum(std::vector<int> v) {
		return std::accumulate(v.begin(),v.end(),0);
	}

	static void proc_counts(brinkman *prob,
		std::vector<int> &lux,std::vector<int> &luy,
		std::vector<int> &lpx,std::vector<int> &lpy);

	static void proc_decomp(brinkman *prob,int &mp,int &np);
};

/**
 * struct containing problem information required
 * for the setup of a Stokes-type FEM problem
 */
struct brinkman {

	/** problem type */
	ptype type;
	petsc_solve::method sol_type;

	// logical grid information

	/** polynomial order */
	int po;
	/** number of elements in horiz/vert directions */
	int m,n;
	/** horizontal / vertical periodicity */
	bool hp,vp;
	/** vel and pressure grids */
	grid gr_u,gr_p;

	// finite element node info

	/** vel and pressure master elements */
	master_quad mast_u,mast_p;
	/** number of velocity nodes in horiz/vert directions */
	int ump,unp;
	/** local,global number of velocity nodes */
	int lu,gu;
	/** number of pressure nodes in horiz/vert directions */
	int pmp,pnp;
	/** local,global number of pressure nodes */
	int lp,gp;

	// problem parameters

	/** pointer to domain struct with access to elem geometry */
	domain *dm;
	/** whether domain is revoled */
	bool rev;
	/** pointers to src_base structs with body force functions */
	tsrc_base<vcomplex> *f;
	/** viscosity */
	double mu;
	/** inverse drag length */
	double alpha;
	/** velocity at infinty */
	vcomplex V_inf;

	/** whether we have data */
	bool data;
	/** velocity and pressure data */
	dcomplex *ux,*uy,*p;

	// constructor/destructor
	brinkman(ptype pt,int po_,int m_,int n_,bool hp_,bool vp_,domain *dm_,
		tsrc_base<vcomplex> *f_,double mu_,double a_,vcomplex V);
	virtual ~brinkman();

	/** whether x-velocity is fixed */
	virtual bool fix_pt_ux(node &n)=0;
	/** whether y-velocity is fixed */
	virtual bool fix_pt_uy(node &n)=0;
	/** fixed x-velocity value */
	virtual dcomplex fix_vl_ux(node &n)=0;
	/** fixed y-velocity value */
	virtual dcomplex fix_vl_uy(node &n)=0;

	// static functions for matrix construction
	static void dm_fe_range(DM dmda,int p,int &eia,int &eib,int &eja,int &ejb);
	static int  node_local_ind(node &n,element &e,DMDALocalInfo &info,int c=0);
	static int  node_natural_ind(node &n,element &e,DMDALocalInfo &info,int c=0);

	// local index retrieval
	void get_u_local_indices(element &e,DM da_vel,int *ind);
	void get_v_local_indices(element &e,DM da_vel,int *ind);
	void get_p_local_indices(element &e,DM da_pre,int *ind);

	// natural index retrieval
	void get_u_natural_indices(element &e,DM da_vel,int *ind);
	void get_v_natural_indices(element &e,DM da_vel,int *ind);
	void get_p_natural_indices(element &e,DM da_pre,int *ind);

	// global index retrieval
	void get_u_global_indices(element &e,DM da_vel,int *ind);
	void get_v_global_indices(element &e,DM da_vel,int *ind);
	void get_p_global_indices(element &e,DM da_pre,int *ind);

	// element calculations
	void get_Auu(element &eu,PetscScalar *m);
	void get_Avv(element &eu,PetscScalar *m);
	void get_Auvp(element &eu,PetscScalar *mx,PetscScalar *my,bool T=false);
	void get_P_mass(element &eu,PetscScalar *m);
	void get_P_stiff(element &eu,PetscScalar *m);
	void get_BU(element &eu,PetscScalar *fx,PetscScalar *fy);

	// matrix block construction
	template<bool x_vel> void add_AUU(element &eu,DM da_vel,Mat Auu);
	void add_Auu(element &eu,DM da_vel,Mat Auu) {add_AUU<true>( eu,da_vel,Auu);}
	void add_Avv(element &eu,DM da_vel,Mat Auu) {add_AUU<false>(eu,da_vel,Auu);}
	void add_Auvp(element &eu,element &ep,DM da_vel,DM da_pre,Mat Aup,Mat Apu);
	void add_P_mass(element &ep,DM da_pre,Mat App);
	void add_P_stiff(element &ep,DM da_pre,Mat App);

	// vector part construction
	template<bool x_vel> void add_Diri_BU(element &eu,DM da_vel,Vec Bu);
	void add_Diri_Bu(element &eu,DM da,Vec Bu) {add_Diri_BU<true>( eu,da,Bu);}
	void add_Diri_Bv(element &eu,DM da,Vec Bu) {add_Diri_BU<false>(eu,da,Bu);}
	void add_Diri_Bp(element &eu,element &ep,DM da_pre,Vec Bp);
	void add_BU(element &eu,DM da_vel,Vec Bu);

	// Dirichlet / fix construction
	void fix_BU(element &eu,DM da_vel,Vec Bu,const dcomplex scale=1);
	void fix_AUU(element &eu,DM da_vel,Mat Auu,const dcomplex scale=1);

	// virtual functions

	/** save problem info in a binary file */
	virtual void save_binary(const char *fn) {bail(-1,"not implemented");}
	/** calculated steady boundary velocity from oscil Brinkman solution */
	virtual vcomplex steady_bvel(node &n) {
		bail(-1,"not implemented"); return vcomplex(0);}
	static brinkman* load_binary(const char *fn);
	static ptype read_type(const char *fn);

	// data management
	void alloc_data();
	void write_data(FILE *fb);
	void read_data(FILE *fb);
	void reset();

	// tools for exact solutions
	virtual vcomplex brink_sol(double x,double y) {
		bail(-1,"not implemented"); return vcomplex(0);};
	static vcomplex brink_exact_sphere(dcomplex V,double a,double alpha,double x,double y);

	static void brink_exact_derivs(dcomplex V, double a,double alpha,double x,double y,
		dcomplex &ux, dcomplex &ux_x,dcomplex &ux_y,dcomplex &ux_xx,dcomplex &ux_xy,dcomplex &ux_yy,
		dcomplex &uy, dcomplex &uy_x,dcomplex &uy_y,dcomplex &uy_xx,dcomplex &uy_xy,dcomplex &uy_yy) {

		double rr = sqrt(x*x + y*y);
		double tt = atan2(y,x);

		dcomplex al = alpha*(1+I)/sqrt(2);

		ux = (V*pow(a,3)*pow(rr,-3))/4. + (3*V*cos(2*tt)*pow(a,3)*pow(rr,-3))/4. \
			  + (3*a*V*pow(al,-2)*pow(rr,-3))/4. + \
			  (9*a*V*cos(2*tt)*pow(al,-2)*pow(rr,-3))/4. + \
			  (3*V*pow(a,2)*pow(al,-1)*pow(rr,-3))/4. + \
			  (9*V*cos(2*tt)*pow(a,2)*pow(al,-1)*pow(rr,-3))/4. - \
			  (3*a*V*pow(al,-2)*exp(a*al - al*rr)*pow(rr,-3))/4. - \
			  (9*a*V*cos(2*tt)*pow(al,-2)*exp(a*al - al*rr)*pow(rr,-3))/4. \
			  - (3*a*V*pow(al,-1)*exp(a*al - al*rr)*pow(rr,-2))/4. - \
			  (9*a*V*cos(2*tt)*pow(al,-1)*exp(a*al - al*rr)*pow(rr,-2))/4. \
			  + (3*a*V*exp(a*al - al*rr)*pow(rr,-1))/4. - \
			  (3*a*V*cos(2*tt)*exp(a*al - al*rr)*pow(rr,-1))/4.;

		ux_x = (3*V*cos(tt)*pow(a,3)*pow(rr,-4))/4. - \
			   (15*V*cos(tt)*cos(2*tt)*pow(a,3)*pow(rr,-4))/4. + \
			   (9*a*V*cos(tt)*pow(al,-2)*pow(rr,-4))/4. - \
			   (45*a*V*cos(tt)*cos(2*tt)*pow(al,-2)*pow(rr,-4))/4. + \
			   (9*V*cos(tt)*pow(a,2)*pow(al,-1)*pow(rr,-4))/4. - \
			   (45*V*cos(tt)*cos(2*tt)*pow(a,2)*pow(al,-1)*pow(rr,-4))/4. - \
			   (9*a*V*cos(tt)*pow(al,-2)*exp(a*al - al*rr)*pow(rr,-4))/4. + \
			   (45*a*V*cos(tt)*cos(2*tt)*pow(al,-2)*exp(a*al - \
				al*rr)*pow(rr,-4))/4. - (9*a*V*cos(tt)*pow(al,-1)*exp(a*al - \
				al*rr)*pow(rr,-3))/4. + \
				(45*a*V*cos(tt)*cos(2*tt)*pow(al,-1)*exp(a*al - \
				al*rr)*pow(rr,-3))/4. - (3*a*V*cos(tt)*exp(a*al - \
				al*rr)*pow(rr,-2))/2. + (9*a*V*cos(tt)*cos(2*tt)*exp(a*al - \
				al*rr)*pow(rr,-2))/2. - (3*a*al*V*cos(tt)*exp(a*al - \
				al*rr)*pow(rr,-1))/4. + (3*a*al*V*cos(tt)*cos(2*tt)*exp(a*al \
				- al*rr)*pow(rr,-1))/4.;

		ux_y = (-9*V*pow(a,3)*pow(rr,-4)*sin(tt))/4. - \
		(15*V*cos(2*tt)*pow(a,3)*pow(rr,-4)*sin(tt))/4. - \
		(27*a*V*pow(al,-2)*pow(rr,-4)*sin(tt))/4. - \
		(45*a*V*cos(2*tt)*pow(al,-2)*pow(rr,-4)*sin(tt))/4. - \
		(27*V*pow(a,2)*pow(al,-1)*pow(rr,-4)*sin(tt))/4. - \
		(45*V*cos(2*tt)*pow(a,2)*pow(al,-1)*pow(rr,-4)*sin(tt))/4. + \
		(27*a*V*pow(al,-2)*exp(a*al - al*rr)*pow(rr,-4)*sin(tt))/4. \
		+ (45*a*V*cos(2*tt)*pow(al,-2)*exp(a*al - \
		al*rr)*pow(rr,-4)*sin(tt))/4. + (27*a*V*pow(al,-1)*exp(a*al \
		- al*rr)*pow(rr,-3)*sin(tt))/4. + \
		(45*a*V*cos(2*tt)*pow(al,-1)*exp(a*al - \
		al*rr)*pow(rr,-3)*sin(tt))/4. + (3*a*V*exp(a*al - \
		al*rr)*pow(rr,-2)*sin(tt))/2. + (9*a*V*cos(2*tt)*exp(a*al - \
		al*rr)*pow(rr,-2)*sin(tt))/2. - (3*a*al*V*exp(a*al - \
		al*rr)*pow(rr,-1)*sin(tt))/4. + (3*a*al*V*cos(2*tt)*exp(a*al \
		- al*rr)*pow(rr,-1)*sin(tt))/4.;

		ux_xx = (27*V*pow(a,3)*pow(rr,-5))/16. + \
		(15*V*cos(2*tt)*pow(a,3)*pow(rr,-5))/4. + \
		(105*V*cos(4*tt)*pow(a,3)*pow(rr,-5))/16. + \
		(81*a*V*pow(al,-2)*pow(rr,-5))/16. + \
		(45*a*V*cos(2*tt)*pow(al,-2)*pow(rr,-5))/4. + \
		(315*a*V*cos(4*tt)*pow(al,-2)*pow(rr,-5))/16. + \
		(81*V*pow(a,2)*pow(al,-1)*pow(rr,-5))/16. + \
		(45*V*cos(2*tt)*pow(a,2)*pow(al,-1)*pow(rr,-5))/4. + \
		(315*V*cos(4*tt)*pow(a,2)*pow(al,-1)*pow(rr,-5))/16. - \
		(81*a*V*pow(al,-2)*exp(a*al - al*rr)*pow(rr,-5))/16. - \
		(45*a*V*cos(2*tt)*pow(al,-2)*exp(a*al - \
				al*rr)*pow(rr,-5))/4. - \
		(315*a*V*cos(4*tt)*pow(al,-2)*exp(a*al - \
				al*rr)*pow(rr,-5))/16. - (81*a*V*pow(al,-1)*exp(a*al - \
				al*rr)*pow(rr,-4))/16. - \
				(45*a*V*cos(2*tt)*pow(al,-1)*exp(a*al - \
				al*rr)*pow(rr,-4))/4. - \
				(315*a*V*cos(4*tt)*pow(al,-1)*exp(a*al - \
				al*rr)*pow(rr,-4))/16. - (33*a*V*exp(a*al - \
				al*rr)*pow(rr,-3))/16. - (9*a*V*cos(2*tt)*exp(a*al - \
				al*rr)*pow(rr,-3))/2. - (135*a*V*cos(4*tt)*exp(a*al - \
				al*rr)*pow(rr,-3))/16. - (3*a*al*V*exp(a*al - \
				al*rr)*pow(rr,-2))/8. - (3*a*al*V*cos(2*tt)*exp(a*al - \
				al*rr)*pow(rr,-2))/4. - (15*a*al*V*cos(4*tt)*exp(a*al - \
				al*rr)*pow(rr,-2))/8. + (3*a*V*pow(al,2)*exp(a*al - \
				al*rr)*pow(rr,-1))/16. - (3*a*V*cos(4*tt)*pow(al,2)*exp(a*al \
				- al*rr)*pow(rr,-1))/16.;

		ux_xy = (15*V*pow(a,3)*pow(rr,-5)*sin(2*tt))/8. + \
				(105*V*cos(2*tt)*pow(a,3)*pow(rr,-5)*sin(2*tt))/8. + \
				(45*a*V*pow(al,-2)*pow(rr,-5)*sin(2*tt))/8. + \
				(315*a*V*cos(2*tt)*pow(al,-2)*pow(rr,-5)*sin(2*tt))/8. + \
				(45*V*pow(a,2)*pow(al,-1)*pow(rr,-5)*sin(2*tt))/8. + \
				(315*V*cos(2*tt)*pow(a,2)*pow(al,-1)*pow(rr,-5)*sin(2*tt))/8. - \
				(45*a*V*pow(al,-2)*exp(a*al - \
										  al*rr)*pow(rr,-5)*sin(2*tt))/8. - \
				(315*a*V*cos(2*tt)*pow(al,-2)*exp(a*al - \
													 al*rr)*pow(rr,-5)*sin(2*tt))/8. - \
				(45*a*V*pow(al,-1)*exp(a*al - \
										  al*rr)*pow(rr,-4)*sin(2*tt))/8. - \
				(315*a*V*cos(2*tt)*pow(al,-1)*exp(a*al - \
					al*rr)*pow(rr,-4)*sin(2*tt))/8. - (9*a*V*exp(a*al - \
					al*rr)*pow(rr,-3)*sin(2*tt))/8. - (135*a*V*cos(2*tt)*exp(a*al - \
					al*rr)*pow(rr,-3)*sin(2*tt))/8. + (3*a*al*V*exp(a*al - \
					al*rr)*pow(rr,-2)*sin(2*tt))/4. - \
				(15*a*al*V*cos(2*tt)*exp(a*al - \
				al*rr)*pow(rr,-2)*sin(2*tt))/4. + (3*a*V*pow(al,2)*exp(a*al \
				- al*rr)*pow(rr,-1)*sin(2*tt))/8. - \
				(3*a*V*cos(2*tt)*pow(al,2)*exp(a*al - \
				al*rr)*pow(rr,-1)*sin(2*tt))/8.;

		ux_yy = (9*V*pow(a,3)*pow(rr,-5))/16. - \
				(105*V*cos(4*tt)*pow(a,3)*pow(rr,-5))/16. + \
				(27*a*V*pow(al,-2)*pow(rr,-5))/16. - \
				(315*a*V*cos(4*tt)*pow(al,-2)*pow(rr,-5))/16. + \
				(27*V*pow(a,2)*pow(al,-1)*pow(rr,-5))/16. - \
				(315*V*cos(4*tt)*pow(a,2)*pow(al,-1)*pow(rr,-5))/16. - \
				(27*a*V*pow(al,-2)*exp(a*al - al*rr)*pow(rr,-5))/16. + \
				(315*a*V*cos(4*tt)*pow(al,-2)*exp(a*al - \
				al*rr)*pow(rr,-5))/16. - (27*a*V*pow(al,-1)*exp(a*al - \
				al*rr)*pow(rr,-4))/16. + \
				(315*a*V*cos(4*tt)*pow(al,-1)*exp(a*al - \
				al*rr)*pow(rr,-4))/16. - (3*a*V*exp(a*al - \
				al*rr)*pow(rr,-3))/16. - (9*a*V*cos(2*tt)*exp(a*al - \
				al*rr)*pow(rr,-3))/4. + (135*a*V*cos(4*tt)*exp(a*al - \
				al*rr)*pow(rr,-3))/16. + (3*a*al*V*exp(a*al - \
				al*rr)*pow(rr,-2))/8. - (9*a*al*V*cos(2*tt)*exp(a*al - \
				al*rr)*pow(rr,-2))/4. + (15*a*al*V*cos(4*tt)*exp(a*al - \
				al*rr)*pow(rr,-2))/8. + (9*a*V*pow(al,2)*exp(a*al - \
				al*rr)*pow(rr,-1))/16. - (3*a*V*cos(2*tt)*pow(al,2)*exp(a*al \
				- al*rr)*pow(rr,-1))/4. + \
				(3*a*V*cos(4*tt)*pow(al,2)*exp(a*al - al*rr)*pow(rr,-1))/16.;

		uy = (3*V*pow(a,3)*pow(rr,-3)*sin(2*tt))/4. + \
			 (9*a*V*pow(al,-2)*pow(rr,-3)*sin(2*tt))/4. + \
			 (9*V*pow(a,2)*pow(al,-1)*pow(rr,-3)*sin(2*tt))/4. - \
			 (9*a*V*pow(al,-2)*exp(a*al - al*rr)*pow(rr,-3)*sin(2*tt))/4. \
			 - (9*a*V*pow(al,-1)*exp(a*al - \
						 al*rr)*pow(rr,-2)*sin(2*tt))/4. - (3*a*V*exp(a*al - \
							 al*rr)*pow(rr,-1)*sin(2*tt))/4.;

		uy_x = (-9*V*pow(a,3)*pow(rr,-4)*sin(tt))/4. - \
			   (15*V*cos(2*tt)*pow(a,3)*pow(rr,-4)*sin(tt))/4. - \
			   (27*a*V*pow(al,-2)*pow(rr,-4)*sin(tt))/4. - \
			   (45*a*V*cos(2*tt)*pow(al,-2)*pow(rr,-4)*sin(tt))/4. - \
			   (27*V*pow(a,2)*pow(al,-1)*pow(rr,-4)*sin(tt))/4. - \
			   (45*V*cos(2*tt)*pow(a,2)*pow(al,-1)*pow(rr,-4)*sin(tt))/4. + \
			   (27*a*V*pow(al,-2)*exp(a*al - al*rr)*pow(rr,-4)*sin(tt))/4. \
			   + (45*a*V*cos(2*tt)*pow(al,-2)*exp(a*al - \
				al*rr)*pow(rr,-4)*sin(tt))/4. + (27*a*V*pow(al,-1)*exp(a*al \
				- al*rr)*pow(rr,-3)*sin(tt))/4. + \
				(45*a*V*cos(2*tt)*pow(al,-1)*exp(a*al - \
				al*rr)*pow(rr,-3)*sin(tt))/4. + 3*a*V*exp(a*al - \
				al*rr)*pow(rr,-2)*sin(tt) + (9*a*V*cos(2*tt)*exp(a*al - \
				al*rr)*pow(rr,-2)*sin(tt))/2. + (3*a*al*V*exp(a*al - \
				al*rr)*pow(rr,-1)*sin(tt))/4. + (3*a*al*V*cos(2*tt)*exp(a*al \
				- al*rr)*pow(rr,-1)*sin(tt))/4.;

		uy_y = (-9*V*cos(tt)*pow(a,3)*pow(rr,-4))/4. + \
				(15*V*cos(tt)*cos(2*tt)*pow(a,3)*pow(rr,-4))/4. - \
				(27*a*V*cos(tt)*pow(al,-2)*pow(rr,-4))/4. + \
				(45*a*V*cos(tt)*cos(2*tt)*pow(al,-2)*pow(rr,-4))/4. - \
				(27*V*cos(tt)*pow(a,2)*pow(al,-1)*pow(rr,-4))/4. + \
				(45*V*cos(tt)*cos(2*tt)*pow(a,2)*pow(al,-1)*pow(rr,-4))/4. + \
				(27*a*V*cos(tt)*pow(al,-2)*exp(a*al - al*rr)*pow(rr,-4))/4. \
				- (45*a*V*cos(tt)*cos(2*tt)*pow(al,-2)*exp(a*al - \
							al*rr)*pow(rr,-4))/4. + (27*a*V*cos(tt)*pow(al,-1)*exp(a*al \
								- al*rr)*pow(rr,-3))/4. - \
							(45*a*V*cos(tt)*cos(2*tt)*pow(al,-1)*exp(a*al - \
						al*rr)*pow(rr,-3))/4. + 3*a*V*cos(tt)*exp(a*al - \
						al*rr)*pow(rr,-2) - (9*a*V*cos(tt)*cos(2*tt)*exp(a*al - \
						al*rr)*pow(rr,-2))/2. + (3*a*al*V*cos(tt)*exp(a*al - \
						al*rr)*pow(rr,-1))/4. - (3*a*al*V*cos(tt)*cos(2*tt)*exp(a*al \
						- al*rr)*pow(rr,-1))/4.;

		uy_xx = (15*V*pow(a,3)*pow(rr,-5)*sin(2*tt))/8. + \
				(105*V*cos(2*tt)*pow(a,3)*pow(rr,-5)*sin(2*tt))/8. + \
				(45*a*V*pow(al,-2)*pow(rr,-5)*sin(2*tt))/8. + \
				(315*a*V*cos(2*tt)*pow(al,-2)*pow(rr,-5)*sin(2*tt))/8. + \
				(45*V*pow(a,2)*pow(al,-1)*pow(rr,-5)*sin(2*tt))/8. + \
				(315*V*cos(2*tt)*pow(a,2)*pow(al,-1)*pow(rr,-5)*sin(2*tt))/8. - \
				(45*a*V*pow(al,-2)*exp(a*al - \
										  al*rr)*pow(rr,-5)*sin(2*tt))/8. - \
				(315*a*V*cos(2*tt)*pow(al,-2)*exp(a*al - \
													 al*rr)*pow(rr,-5)*sin(2*tt))/8. - \
				(45*a*V*pow(al,-1)*exp(a*al - \
										  al*rr)*pow(rr,-4)*sin(2*tt))/8. - \
				(315*a*V*cos(2*tt)*pow(al,-1)*exp(a*al - \
				al*rr)*pow(rr,-4)*sin(2*tt))/8. - (27*a*V*exp(a*al - \
				al*rr)*pow(rr,-3)*sin(2*tt))/8. - (135*a*V*cos(2*tt)*exp(a*al - \
				al*rr)*pow(rr,-3)*sin(2*tt))/8. - (3*a*al*V*exp(a*al - \
				al*rr)*pow(rr,-2)*sin(2*tt))/2. - \
				(15*a*al*V*cos(2*tt)*exp(a*al - \
				al*rr)*pow(rr,-2)*sin(2*tt))/4. - (3*a*V*pow(al,2)*exp(a*al \
				- al*rr)*pow(rr,-1)*sin(2*tt))/8. - \
				(3*a*V*cos(2*tt)*pow(al,2)*exp(a*al - \
				al*rr)*pow(rr,-1)*sin(2*tt))/8.;

		uy_xy = (9*V*pow(a,3)*pow(rr,-5))/16. - \
				 (105*V*cos(4*tt)*pow(a,3)*pow(rr,-5))/16. + \
				 (27*a*V*pow(al,-2)*pow(rr,-5))/16. - \
				 (315*a*V*cos(4*tt)*pow(al,-2)*pow(rr,-5))/16. + \
				 (27*V*pow(a,2)*pow(al,-1)*pow(rr,-5))/16. - \
				 (315*V*cos(4*tt)*pow(a,2)*pow(al,-1)*pow(rr,-5))/16. - \
				 (27*a*V*pow(al,-2)*exp(a*al - al*rr)*pow(rr,-5))/16. + \
				 (315*a*V*cos(4*tt)*pow(al,-2)*exp(a*al - \
				al*rr)*pow(rr,-5))/16. - (27*a*V*pow(al,-1)*exp(a*al - \
				al*rr)*pow(rr,-4))/16. + \
				(315*a*V*cos(4*tt)*pow(al,-1)*exp(a*al - \
				al*rr)*pow(rr,-4))/16. - (15*a*V*exp(a*al - \
				al*rr)*pow(rr,-3))/16. + (135*a*V*cos(4*tt)*exp(a*al - \
				al*rr)*pow(rr,-3))/16. - (3*a*al*V*exp(a*al - \
				al*rr)*pow(rr,-2))/8. + (15*a*al*V*cos(4*tt)*exp(a*al - \
				al*rr)*pow(rr,-2))/8. - (3*a*V*pow(al,2)*exp(a*al - \
				al*rr)*pow(rr,-1))/16. + (3*a*V*cos(4*tt)*pow(al,2)*exp(a*al \
				- al*rr)*pow(rr,-1))/16.;

		uy_yy = (15*V*pow(a,3)*pow(rr,-5)*sin(2*tt))/8. - \
				(105*V*cos(2*tt)*pow(a,3)*pow(rr,-5)*sin(2*tt))/8. + \
				(45*a*V*pow(al,-2)*pow(rr,-5)*sin(2*tt))/8. - \
				(315*a*V*cos(2*tt)*pow(al,-2)*pow(rr,-5)*sin(2*tt))/8. + \
				(45*V*pow(a,2)*pow(al,-1)*pow(rr,-5)*sin(2*tt))/8. - \
				(315*V*cos(2*tt)*pow(a,2)*pow(al,-1)*pow(rr,-5)*sin(2*tt))/8. - \
				(45*a*V*pow(al,-2)*exp(a*al - \
										  al*rr)*pow(rr,-5)*sin(2*tt))/8. + \
				(315*a*V*cos(2*tt)*pow(al,-2)*exp(a*al - \
													 al*rr)*pow(rr,-5)*sin(2*tt))/8. - \
				(45*a*V*pow(al,-1)*exp(a*al - \
										  al*rr)*pow(rr,-4)*sin(2*tt))/8. + \
				(315*a*V*cos(2*tt)*pow(al,-1)*exp(a*al - \
				al*rr)*pow(rr,-4)*sin(2*tt))/8. - (27*a*V*exp(a*al - \
				al*rr)*pow(rr,-3)*sin(2*tt))/8. + (135*a*V*cos(2*tt)*exp(a*al - \
				al*rr)*pow(rr,-3)*sin(2*tt))/8. - (3*a*al*V*exp(a*al - \
				al*rr)*pow(rr,-2)*sin(2*tt))/2. + \
				(15*a*al*V*cos(2*tt)*exp(a*al - \
				al*rr)*pow(rr,-2)*sin(2*tt))/4. - (3*a*V*pow(al,2)*exp(a*al \
				- al*rr)*pow(rr,-1)*sin(2*tt))/8. + \
				(3*a*V*cos(2*tt)*pow(al,2)*exp(a*al - \
				al*rr)*pow(rr,-1)*sin(2*tt))/8.;

		/*
		ux_x = (-3*a*V*cos(tt)*pow(alpha,-2)*exp(-(alpha*rr))*pow(rr,-4)*(-((3 + 3*a*alpha + pow(a,2)*pow(alpha,2))*exp(alpha*rr)) + exp(a*alpha)*(3 + 3*alpha*rr + 2*pow(alpha,2)*pow(rr,2) + pow(alpha,3)*pow(rr,3)) - cos(2*tt)*(-5*(3 + 3*a*alpha + pow(a,2)*pow(alpha,2))*exp(alpha*rr) + exp(a*alpha)*(15 + 15*alpha*rr + 6*pow(alpha,2)*pow(rr,2) + pow(alpha,3)*pow(rr,3)))))/4.;
		ux_y = (-3*a*V*pow(alpha,-2)*exp(-(alpha*rr))*pow(rr,-4)*(3*(3 + 3*a*alpha + pow(a,2)*pow(alpha,2))*exp(alpha*rr) + exp(a*alpha)*(-9 - 9*alpha*rr - 2*pow(alpha,2)*pow(rr,2) + pow(alpha,3)*pow(rr,3)) - cos(2*tt)*(-5*(3 + 3*a*alpha + pow(a,2)*pow(alpha,2))*exp(alpha*rr) + exp(a*alpha)*(15 + 15*alpha*rr + 6*pow(alpha,2)*pow(rr,2) + pow(alpha,3)*pow(rr,3))))*sin(tt))/4.;
		ux_xx = (3*a*V*pow(alpha,-2)*exp(-(alpha*rr))*pow(rr,-5)*(-27*exp(a*alpha) - 27*alpha*rr*exp(a*alpha) + 27*exp(alpha*rr) + 27*a*alpha*exp(alpha*rr) + 9*pow(a,2)*pow(alpha,2)*exp(alpha*rr) - 11*pow(alpha,2)*exp(a*alpha)*pow(rr,2) - 2*pow(alpha,3)*exp(a*alpha)*pow(rr,3) - 4*cos(2*tt)*(-5*(3 + 3*a*alpha + pow(a,2)*pow(alpha,2))*exp(alpha*rr) + exp(a*alpha)*(15 + 15*alpha*rr + 6*pow(alpha,2)*pow(rr,2) + pow(alpha,3)*pow(rr,3))) + pow(alpha,4)*exp(a*alpha)*pow(rr,4) - cos(4*tt)*(-35*(3 + 3*a*alpha + pow(a,2)*pow(alpha,2))*exp(alpha*rr) + exp(a*alpha)*(105 + 105*alpha*rr + 45*pow(alpha,2)*pow(rr,2) + 10*pow(alpha,3)*pow(rr,3) + pow(alpha,4)*pow(rr,4)))))/16.;
		ux_xy = (3*a*V*pow(alpha,-2)*exp(-(alpha*rr))*pow(rr,-5)*(5*(3 + 3*a*alpha + pow(a,2)*pow(alpha,2))*exp(alpha*rr) + exp(a*alpha)*(-15 - 15*alpha*rr - 3*pow(alpha,2)*pow(rr,2) + 2*pow(alpha,3)*pow(rr,3) + pow(alpha,4)*pow(rr,4)) - cos(2*tt)*(-35*(3 + 3*a*alpha + pow(a,2)*pow(alpha,2))*exp(alpha*rr) + exp(a*alpha)*(105 + 105*alpha*rr + 45*pow(alpha,2)*pow(rr,2) + 10*pow(alpha,3)*pow(rr,3) + pow(alpha,4)*pow(rr,4))))*sin(2*tt))/8.;
		ux_yy = (-3*a*V*pow(alpha,-2)*exp(-(alpha*rr))*pow(rr,-5)*(9*exp(a*alpha) + 9*alpha*rr*exp(a*alpha) - 9*exp(alpha*rr) - 9*a*alpha*exp(alpha*rr) - 3*pow(a,2)*pow(alpha,2)*exp(alpha*rr) + pow(alpha,2)*exp(a*alpha)*pow(rr,2) + 4*cos(2*tt)*pow(alpha,2)*exp(a*alpha)*pow(rr,2)*(3 + 3*alpha*rr + pow(alpha,2)*pow(rr,2)) - 2*pow(alpha,3)*exp(a*alpha)*pow(rr,3) - 3*pow(alpha,4)*exp(a*alpha)*pow(rr,4) - cos(4*tt)*(-35*(3 + 3*a*alpha + pow(a,2)*pow(alpha,2))*exp(alpha*rr) + exp(a*alpha)*(105 + 105*alpha*rr + 45*pow(alpha,2)*pow(rr,2) + 10*pow(alpha,3)*pow(rr,3) + pow(alpha,4)*pow(rr,4)))))/16.;

		uy_x =(-3*a*V*pow(alpha,-2)*exp(-(alpha*rr))*pow(rr,-4)*(3*(3 + 3*a*alpha + pow(a,2)*pow(alpha,2))*exp(alpha*rr) - exp(a*alpha)*(9 + 9*alpha*rr + 4*pow(alpha,2)*pow(rr,2) + pow(alpha,3)*pow(rr,3)) - cos(2*tt)*(-5*(3 + 3*a*alpha + pow(a,2)*pow(alpha,2))*exp(alpha*rr) + exp(a*alpha)*(15 + 15*alpha*rr + 6*pow(alpha,2)*pow(rr,2) + pow(alpha,3)*pow(rr,3))))*sin(tt))/4.;
		uy_y = (3*a*V*cos(tt)*pow(alpha,-2)*exp(-(alpha*rr))*pow(rr,-4)*(-3*(3 + 3*a*alpha + pow(a,2)*pow(alpha,2))*exp(alpha*rr) + exp(a*alpha)*(9 + 9*alpha*rr + 4*pow(alpha,2)*pow(rr,2) + pow(alpha,3)*pow(rr,3)) - cos(2*tt)*(-5*(3 + 3*a*alpha + pow(a,2)*pow(alpha,2))*exp(alpha*rr) + exp(a*alpha)*(15 + 15*alpha*rr + 6*pow(alpha,2)*pow(rr,2) + pow(alpha,3)*pow(rr,3)))))/4.;
		uy_xx = (3*a*V*pow(alpha,-2)*exp(-(alpha*rr))*pow(rr,-5)*(5*(3 + 3*a*alpha + pow(a,2)*pow(alpha,2))*exp(alpha*rr) - exp(a*alpha)*(15 + 15*alpha*rr + 9*pow(alpha,2)*pow(rr,2) + 4*pow(alpha,3)*pow(rr,3) + pow(alpha,4)*pow(rr,4)) - cos(2*tt)*(-35*(3 + 3*a*alpha + pow(a,2)*pow(alpha,2))*exp(alpha*rr) + exp(a*alpha)*(105 + 105*alpha*rr + 45*pow(alpha,2)*pow(rr,2) + 10*pow(alpha,3)*pow(rr,3) + pow(alpha,4)*pow(rr,4))))*sin(2*tt))/8.;
		uy_xy = (-3*a*V*pow(alpha,-2)*exp(-(alpha*rr))*pow(rr,-5)*(-3*(3 + 3*a*alpha + pow(a,2)*pow(alpha,2))*exp(alpha*rr) + exp(a*alpha)*(9 + 9*alpha*rr + 5*pow(alpha,2)*pow(rr,2) + 2*pow(alpha,3)*pow(rr,3) + pow(alpha,4)*pow(rr,4)) - cos(4*tt)*(-35*(3 + 3*a*alpha + pow(a,2)*pow(alpha,2))*exp(alpha*rr) + exp(a*alpha)*(105 + 105*alpha*rr + 45*pow(alpha,2)*pow(rr,2) + 10*pow(alpha,3)*pow(rr,3) + pow(alpha,4)*pow(rr,4)))))/16.;
		uy_yy = (-3*a*V*pow(alpha,-2)*exp(-(alpha*rr))*pow(rr,-5)*(-5*(3 + 3*a*alpha + pow(a,2)*pow(alpha,2))*exp(alpha*rr) + exp(a*alpha)*(15 + 15*alpha*rr + 9*pow(alpha,2)*pow(rr,2) + 4*pow(alpha,3)*pow(rr,3) + pow(alpha,4)*pow(rr,4)) - cos(2*tt)*(-35*(3 + 3*a*alpha + pow(a,2)*pow(alpha,2))*exp(alpha*rr) + exp(a*alpha)*(105 + 105*alpha*rr + 45*pow(alpha,2)*pow(rr,2) + 10*pow(alpha,3)*pow(rr,3) + pow(alpha,4)*pow(rr,4))))*sin(2*tt))/8.;

		*/

	}

	// tools for coarsening
	brinkman* coarsen(DM comp);

	enum class dmgr {sys,vel,pre};
	static dmgr dm_type(DM dm) {
		const char *op;
		DMGetOptionsPrefix(dm,&op);

		// get last 4 chars
		char pref[5];
		int len = strlen(op);
		if (len < 4) {
			bail(-1,"options prefix has too few chars (%d)\n",len);
		}

		char *pp=pref;
		for(int i=len-4;i<=len;i++) (*pp++)=op[i];

		if      (strcmp(pref,"vel_")==0) return dmgr::vel;
		else if (strcmp(pref,"pre_")==0) return dmgr::pre;
		else if (strcmp(pref,"sys_")==0) return dmgr::sys;
		else       bail(-1,"bad option prefix \"%s\"\n",pref);

		return dmgr::sys;
	}
	void coarsen_factors_vel(int Mu,int Nu,int &cx,int &cy);
	void coarsen_factors_pre(int Mp,int Np,int &cx,int &cy);
	void coarsen_factors(int Mu,int Nu,int Mp,int Np,int &cx,int &cy);
	virtual brinkman* coarsen(int cx,int cy) {
		bail(-1,"not implemented"); return NULL;};
	
	// reading data from petsc
	void read(petsc_solve &solver,Vec X=NULL);

	// solve with petsc
	void solve() {
		petsc_solve solver(this,sol_type);
		solver.solve();
	}

	// solution accessors
	
	// interpolated solution
	template<dcomplex* brinkman::*arr>
	void interp_sol(double x,double y,dcomplex &v,dcomplex &v_x,dcomplex &v_y);

	// interpolated solution and derivatives
	template<dcomplex* brinkman::*arr>
	dcomplex interp_sol(double x,double y);

	// fields
	vcomplex get_u(double x,double y){
		return vcomplex(get_ux(x,y),get_uy(x,y));
	}
	dcomplex get_ux(double x,double y) {return interp_sol<&brinkman::ux>(x,y);}
	dcomplex get_uy(double x,double y) {return interp_sol<&brinkman::uy>(x,y);}
	dcomplex get_p( double x,double y) {return interp_sol<&brinkman::p >(x,y);}

	void get_u_derivs(double x,double y,dcomplex &ux,dcomplex &ux_x,dcomplex &ux_y,
			dcomplex &ux_xx,dcomplex &ux_xy,dcomplex &ux_yy) {

		// reference to master element
		master_quad &mast = mast_u;

		// reference to grid
		grid &gr = gr_u;

		// number of nodes in an element
		int ln = lu;
		dcomplex* val = new dcomplex[ln];

		// requested solution array
		dcomplex* sol = this->ux;

		// grab element index and local coords
		int ei,ej;
		double xi,eta;


		// grab element if we got one in the grid
		if (dm->inv_coords(x,y,ei,ej,xi,eta)) {

			element e = gr(ei,ej);

			trafo_base *ct = dm->alloc_element_geometry(e);


			// fill vals array
			for(node n:e) val[n.l]=sol[n.g];

			// evaluate result
			mast.eval<dcomplex>(*ct,xi,eta,val,
				ux,ux_x,ux_y,ux_xx,ux_xy,ux_yy
					);

			dm->free_element_geometry(ct);

		// otherwise, a nan
		} else {
			ux=ux_x=ux_y=ux_xx=ux_xy=ux_yy = std::nan("out of grid")*(1+I);
		}

		// cleanup and return
		delete[] val;
	}

	void get_v_derivs(double x,double y,dcomplex &uy,dcomplex &uy_x,dcomplex &uy_y,
			dcomplex &uy_xx,dcomplex &uy_xy,dcomplex &uy_yy) {

		// reference to master element
		master_quad &mast = mast_u;

		// reference to grid
		grid &gr = gr_u;

		// number of nodes in an element
		int ln = lu;
		dcomplex* val = new dcomplex[ln];

		// requested solution array
		dcomplex* sol = this->uy;

		// grab element index and local coords
		int ei,ej;
		double xi,eta;


		// grab element if we got one in the grid
		if (dm->inv_coords(x,y,ei,ej,xi,eta)) {

			element e = gr(ei,ej);

			trafo_base *ct = dm->alloc_element_geometry(e);


			// fill vals array
			for(node n:e) val[n.l]=sol[n.g];

			// evaluate result
			mast.eval<dcomplex>(*ct,xi,eta,val,
				uy,uy_x,uy_y,uy_xx,uy_xy,uy_yy
					);

			dm->free_element_geometry(ct);

		// otherwise, a nan
		} else {
			uy=uy_x=uy_y=uy_xx=uy_xy=uy_yy = std::nan("out of grid")*(1+I);
		}

		// cleanup and return
		delete[] val;
	}

	void get_p_derivs(double x,double y,dcomplex &p,dcomplex &p_x,dcomplex &p_y,
			dcomplex &p_xx,dcomplex &p_xy,dcomplex &p_yy) {

		// reference to master element
		master_quad &mast = mast_p;

		// reference to grid
		grid &gr = gr_p;

		// number of nodes in an element
		int ln = lp;
		dcomplex* val = new dcomplex[ln];

		// requested solution array
		dcomplex* sol = this->p;

		// grab element index and local coords
		int ei,ej;
		double xi,eta;


		// grab element if we got one in the grid
		if (dm->inv_coords(x,y,ei,ej,xi,eta)) {

			element e = gr(ei,ej);

			trafo_base *ct = dm->alloc_element_geometry(e);


			// fill vals array
			for(node n:e) val[n.l]=sol[n.g];

			// evaluate result
			mast.eval<dcomplex>(*ct,xi,eta,val,
				p,p_x,p_y,p_xx,p_xy,p_yy
					);

			dm->free_element_geometry(ct);

		// otherwise, a nan
		} else {
			p=p_x=p_y=p_xx=p_xy=p_yy = std::nan("out of grid")*(1+I);
		}

		// cleanup and return
		delete[] val;
	}

	// deriv operators
	dcomplex get_w(double x,double y); // vorticity
	dcomplex get_d(double x,double y); // divergence

	// "exact" solutions
	dcomplex get_u_brink(double x,double y) {return brink_sol(x,y).x;}
	dcomplex get_v_brink(double x,double y) {return brink_sol(x,y).y;}

	// steady streaming part of advective term
	vcomplex u_grad_u(double x,double y);
	vcomplex u_grad_u(element &e,node &n);

	// gnuplot output

	// generalized
	template<int t,bool real> void gp_quad(const char* fn);
	template<int t,bool real> void gp_line(const char* fn);

	/** real x-velocity quads */
	void ux_r_quad(const char* fn) {gp_quad<2,true>(fn);}
	/** real y-velocity quads */
	void uy_r_quad(const char* fn) {gp_quad<1,true>(fn);}
	/** real pressure quads */
	void p_r_quad(const char* fn)  {gp_quad<0,true>(fn);}
	/** real x-velocity lines */
	void ux_r_line(const char* fn) {gp_line<2,true>(fn);}
	/** real y-velocity lines */
	void uy_r_line(const char* fn) {gp_line<1,true>(fn);}
	/** real pressure lines */
	void p_r_line(const char* fn)  {gp_line<0,true>(fn);}
	/** imag x-velocity quads */
	void ux_i_quad(const char* fn) {gp_quad<2,false>(fn);}
	/** imag y-velocity quads */
	void uy_i_quad(const char* fn) {gp_quad<1,false>(fn);}
	/** imag pressure quads */
	void p_i_quad(const char* fn)  {gp_quad<0,false>(fn);}
	/** imag x-velocity lines */
	void ux_i_line(const char* fn) {gp_line<2,false>(fn);}
	/** imag y-velocity lines */
	void uy_i_line(const char* fn) {gp_line<1,false>(fn);}
	/** imag pressure lines */
	void p_i_line(const char* fn)  {gp_line<0,false>(fn);}
};

/**
 * Sample problem type on a simple rectangular grid, useful for testing
 * performance of geometric multigrid
 */
struct driven_cyl: public brinkman {

	/** whether to delete source struct */
	bool del_src;
	/** height of box */
	double h;
	/** length of box */
	double L;

	// base constructor
	driven_cyl(int p,int m,int n,bool hp,bool vp,double h_,double L_,
		tsrc_base<vcomplex> *f,double mu,double alph,vcomplex V_):
		brinkman(ptype::driven_cyl,p,m,n,hp,vp,
			new cartesian_grid(m,n,-L_/2,L_/2,-h_/2,h_/2),f,mu,alph,V_),
			del_src(true),h(h_),L(L_) {}
	driven_cyl(int p,int m,int n,double h,double L,double mu,double alpha,double V):
		driven_cyl(p,m,n,false,false,h,L,new homog(),mu,alpha,vcomplex(V,0)) {}

	// destructor
	~driven_cyl() {if (del_src) delete f;}

	/** whether a x-velocity is fixed (sphere surfaces and inf) */
	bool fix_pt_ux(node &n) {return n.gi==0||n.gi==(ump-1)||n.gj==0||n.gj==(unp-1);}
	/** fixed x-velocity value (1 on spheres, 0 at inf) */
	dcomplex fix_vl_ux(node &n) {return (n.gi==0||n.gi==(ump-1))?0:1;}
	/** whether a y-velocity is fixed (sphere surfaces and inf) */
	bool fix_pt_uy(node &n) {return n.gj==0||n.gj==(unp-1)||n.gi==0||n.gi==(ump-1);}
	/** fixed y-velocity value (zero) */
	dcomplex fix_vl_uy(node &n){return 0;}

	// coarsening
	brinkman *coarsen(int cx,int cy);

	// saving and loading
	void save_binary(const char *fn) {bail(-1,"don't do it\n");}
	static driven_cyl* load_binary(const char *fn) {
		bail(-1,"don't do it\n"); return NULL;}
};

/** struct corresponding to 3D axisymmetric spherical problem */
struct one_sphere: public brinkman {

	/** whether to delete source */
	bool del_src;
	/** radius of sphere */
	double a;
	/** outer radius of domain */
	double R;
	/** const x-,y-velocity of inner sphere (outer is inf)  */
	vcomplex V_i;
	/** varying x-,y-velocity of inner, outer sphere */
	vcomplex *dV_i,*dV_o;

	// base constructor
	one_sphere(int p,int m,int n,
		double a,double R,tsrc_base<vcomplex> *f,
		double mu,double alph,vcomplex Vi,vcomplex Vo);

	// single-velocity constructor
	one_sphere(int p,int m,int n,
		double a,double R,double mu,double al,dcomplex V):
		one_sphere(p,m,n,a,R,new homog(),mu,al,vcomplex(V,0),vcomplex::zero) {}

	// steady / oscillatory constructors
	one_sphere(one_sphere *p,int type=0);

	// destructor
	~one_sphere();

	// boundary velocity due to steady streaming
	vcomplex steady_bvel(node &n);

	// helpers for where we are in the domain
	bool in(node &n) {return n.gi==0;}
	bool out(node &n) {return n.gi==ump-1;}
	bool yaxis(node &n) {return n.gj==0||n.gj==(unp/2);}

	// helpers for vels on spheres
	dcomplex in_vx(int i)  {return   V_i.x + dV_i[i].x;}
	dcomplex out_vx(int i) {return V_inf.x + dV_o[i].x;}
	dcomplex in_vy(int i)  {return   V_i.y + dV_i[i].y;}
	dcomplex out_vy(int i) {return V_inf.y + dV_o[i].y;}

	/** whether a x-velocity is fixed (sphere surfaces and inf) */
	bool fix_pt_ux(node &n) {return in(n) || out(n);}
	/** fixed x-velocity value (1 on spheres, 0 at inf) */
	dcomplex fix_vl_ux(node &n) {return in(n) ? in_vx(n.gj) : out_vx(n.gj);}
	/** whether a y-velocity is fixed (sphere surfaces and inf) */
	bool fix_pt_uy(node &n) {return in(n) || out(n) || yaxis(n);}
	/** fixed y-velocity value (zero) */
	dcomplex fix_vl_uy(node &n){return in(n)?in_vy(n.gj):(out(n)?out_vy(n.gj):0);}

	// save / load binaries
	void save_binary(const char *fn);
	static one_sphere* load_binary(const char *fn);

	// one sphere exact solution
	vcomplex brink_sol(double x,double y);

	// coarsening
	brinkman *coarsen(int cx,int cy);
};


/** struct corresponding to the 3D axisymmetric two-spheres problem */
struct two_spheres: public brinkman {

	/** whether to delete source */
	bool del_src;
	/** radius of left and right spheres */
	double rl,rr;
	/** sphere separation */
	double d;
	/** x-coords of left,right sphere */
	double xl,xr;
	/** xi-coords of left,right sphere */
	double xil,xir;
	/** sphere density (relative to fluid) */
	double rho_l,rho_r;
	/** const x-,y-velocity of left sphere */
	vcomplex V_l;
	/** const x-,y-velocity of right sphere */
	vcomplex V_r;
	/** varying x-,y-velocity of left sphere */
	vcomplex *dV_l;
	/** varying x-,y-velocity of right sphere */
	vcomplex *dV_r;
	/** applied force over both spheres */
	vcomplex F_a;

	// base constructor
	two_spheres(int p,int m,int n,
		double D,double r1,double r2,tsrc_base<vcomplex> *f,double mu,
		double alph,vcomplex V_l,vcomplex V_r,vcomplex V_inf);

	// single vel (actually vel ratio) constructor
	two_spheres(int p,int m,int n,double D,double r1,double r2,
		double mu,double alph,dcomplex Vl,dcomplex Vr):
		two_spheres(p,m,n,D,r1,r2,new homog(),mu,alph,
		vcomplex(Vl,0),vcomplex(Vr,0),vcomplex::zero) {}

	// stokes towing constructor
	two_spheres(int p,int m,int n,double D,double r1,double r2,double mu,
		dcomplex Vtow):two_spheres(p,m,n,D,r1,r2,mu,0,Vtow,Vtow) {}

	// steady stokes constructors
	two_spheres(two_spheres*,dcomplex V=0,int type=0);

	// destructor
	~two_spheres();

	/** set sphere densities */
	void set_densities(double rl,double rr) {rho_l=rl; rho_r=rr;}

	// boundary velocity due to steady streaming
	vcomplex steady_bvel(node &n);

	// geometric helper functions
	bool inf(node &n) {return n.gj==(unp/2) && n.gi==(ump/2);}
	bool left(node &n) {return n.gi==0;}
	bool right(node &n) {return n.gi==(ump-1);}
	bool yaxis(node &n) {return n.gj==0 || n.gj==(unp/2);}

	dcomplex l_vx(int i) {return V_l.x + dV_l[i].x;}
	dcomplex r_vx(int i) {return V_r.x + dV_r[i].x;}
	dcomplex l_vy(int i) {return V_l.y + dV_l[i].y;}
	dcomplex r_vy(int i) {return V_r.y + dV_r[i].y;}

	// helper for setting <Vx-Vr> = 1
	void norm_Vx() {

		// get magnitude of oscillation magnitude
		dcomplex sum(V_r.x-V_l.x);
		//sum *= std::conj(sum);
		//sum = sqrt(std::real(sum));

		// renormalize
		V_r/=sum; V_l/=sum;
	}

	/** whether a x-velocity is fixed (sphere surfaces and inf) */
	bool fix_pt_ux(node &n) {return left(n)||right(n)||inf(n);}
	/** fixed x-velocity value (1 on spheres, 0 at inf) */
	dcomplex fix_vl_ux(node &n) {return inf(n) ? V_inf.x :
		(left(n) ? l_vx(n.gj) : r_vx(n.gj));}
	/** whether a y-velocity is fixed (sphere surfaces and inf) */
	bool fix_pt_uy(node &n) {return left(n)||right(n)||inf(n)||yaxis(n);}
	/** fixed y-velocity value (zero) */
	dcomplex fix_vl_uy(node &n) {return (inf(n) ? V_inf.y : (yaxis(n) ? 0 :
		(left(n) ? l_vy(n.gj) : r_vy(n.gj))));}

	// save/load binary
	void save_binary(const char *fn);
	static two_spheres* load_binary(const char *fn);

	// two spheres exact solution
	vcomplex brink_sol(double x,double y);
	void brink_derivs(double x,double y,
			dcomplex &ux,dcomplex &ux_x,dcomplex &ux_y,dcomplex &ux_xx,dcomplex &ux_xy,dcomplex &ux_yy,
			dcomplex &uy,dcomplex &uy_x,dcomplex &uy_y,dcomplex &uy_xx,dcomplex &uy_xy,dcomplex &uy_yy) {

		// check for x and y nan (at inf)
		if ((std::isnan(x) || std::isinf(x)) && (std::isnan(y) || std::isinf(y))) {
			ux_x = ux_y = ux_xx = ux_xy = ux_yy = uy_x = uy_y
				= uy_xx = uy_xy = uy_yy = 0;
			return;
		}



		dcomplex ux_l,ux_x_l,ux_y_l,ux_xx_l,ux_xy_l,ux_yy_l;
		dcomplex ux_r,ux_x_r,ux_y_r,ux_xx_r,ux_xy_r,ux_yy_r;
		dcomplex uy_l,uy_x_l,uy_y_l,uy_xx_l,uy_xy_l,uy_yy_l;
		dcomplex uy_r,uy_x_r,uy_y_r,uy_xx_r,uy_xy_r,uy_yy_r;

		brinkman::brink_exact_derivs(V_l.x, rl, alpha, x-xl, y,
				ux_l,ux_x_l, ux_y_l, ux_xx_l, ux_xy_l, ux_yy_l,
				uy_l,uy_x_l, uy_y_l, uy_xx_l, uy_xy_l, uy_yy_l);

		brinkman::brink_exact_derivs(V_r.x, rr, alpha, x-xr, y,
				ux_r,ux_x_r, ux_y_r, ux_xx_r, ux_xy_r, ux_yy_r,
				uy_r,uy_x_r, uy_y_r, uy_xx_r, uy_xy_r, uy_yy_r);

		ux = ux_l + ux_r;
		uy = uy_l + uy_r;
		ux_x = ux_x_l + ux_x_r;
		ux_y = ux_y_l + ux_y_r;
		ux_xx = ux_xx_l + ux_xx_r;
		ux_xy = ux_xy_l + ux_xy_r;
		ux_yy = ux_yy_l + ux_yy_r;

		uy_x = uy_x_l + uy_x_r;
		uy_y = uy_y_l + uy_y_r;
		uy_xx = uy_xx_l + uy_xx_r;
		uy_xy = uy_xy_l + uy_xy_r;
		uy_yy = uy_yy_l + uy_yy_r;
	}

	dcomplex avg_vel(bool left);
	// force integration due to oscillating traction
	dcomplex osc_force(bool left);
	// force-free integration functions
	dcomplex drag_force(bool left,bool approx=false);
	dcomplex drag_coeff(bool approx=false);
	/** applied force necessary for accel + drag */
	dcomplex appl_force(bool left) {
		return net_force(left) - drag_force(left);
	}
	/** net force calculated from acceleration */
	dcomplex net_force(bool left) {

		dcomplex V = left?V_l.x:V_r.x;
		double r = left?rl:rr;

		return (I*4*M_PI/3) * r*r*r * V * (left?rho_l:rho_r);

		dcomplex v;
		double sep = xr-xl;
		if (left) {
			v =   V_r.x * pow(rr/sep,3);
		} else {
			v =   V_l.x * pow(rl/sep,3);
		}

		return (I*4*M_PI/3) * r*r*r * V * (left?rho_l:rho_r) +
			(I*2*M_PI/3) * r*r*r * (V-v) * (left?rho_l:rho_r);
	}

	/** initialize boundary velocities to solve for a known solution */
	void init_test_boundary() {

		trafo_base *ct_l, *ct_r;

		double xi_l,et_l, xi_r,et_r;
		double x_l,y_l, x_r,y_r;

		// for each element along sphere surfaces
		for (int gj=0, gi_l=0, gi_r=(ump-1); gj<unp; gj++) {

			// left and right nodes
			node n_l(     (gj%po) * (po+1),  0, gj%po, gi_l + gj*ump, gi_l, gj);
			node n_r(po + (gj%po) * (po+1), po, gj%po, gi_r + gj*ump, gi_r, gj);

			// left and right element nodes
			int ei_l = 0, ei_r = m-1;
			int ej = gj/po;

			element e_l = gr_u(ei_l,ej), e_r = gr_u(ei_r,ej);

			ct_l = dm->alloc_element_geometry(e_l);
			ct_r = dm->alloc_element_geometry(e_r);

			// local coords

//			printf("n_l.l=%d, n_r.l=%d\n",n_l.l,n_r.l);

			mast_u.point(n_l.l, xi_l, et_l);
			mast_u.point(n_r.l, xi_r, et_r);

			// cartesian coords
			ct_l->xi2x(xi_l,et_l,x_l,y_l);
			ct_r->xi2x(xi_r,et_r,x_r,y_r);

			// 
			vcomplex v_l = brink_sol(x_l,y_l);
			vcomplex v_r = brink_sol(x_r,y_r);

			dV_l[gj] = v_l-V_l;
			dV_r[gj] = v_r-V_r;

			dm->free_element_geometry(ct_l);
			dm->free_element_geometry(ct_r);

//			printf("gj = %d, (x_l,y_l)=(%g,%g), (x_r,y_r)=(%g,%g), dV_lX[%d] = %g + i*%g, dV_rX[%d] = %g + i*%g\n",
//					gj,x_l,y_l,x_r,y_r,
//					gj,std::real(dV_l[gj].x),std::imag(dV_l[gj].x),std::real(dV_r[gj].x),std::imag(dV_r[gj].x));
		}
	}

	void calc_boundary_error(FILE *fh,bool left=false) {

		// 1D gauss rule
		gauss_quad_1d q((gr_u.p+1)/2 + 3);

		// coord trafos
		trafo_base *ct_u;

		dcomplex ux_v,ux_x,ux_y,ux_xx,ux_xy,ux_yy;
		dcomplex uy_v,uy_x,uy_y,uy_xx,uy_xy,uy_yy;
		dcomplex ux_k_v,ux_k_x,ux_k_y,ux_k_xx,ux_k_xy,ux_k_yy;
		dcomplex uy_k_v,uy_k_x,uy_k_y,uy_k_xx,uy_k_xy,uy_k_yy;

		// solution datapoints
		dcomplex *usol = new dcomplex[lu],
			*vsol = new dcomplex[lu];

		// Jacobean info / normal vector
		double det;
		double x,y;
		double xi_x,et_x,xi_y,et_y,x_et,y_et;
		double nx,ny,n2,nx_x,ny_x;

		//mesg("%s\n",left?"left!!!":"right!!!");

		double err_v=0;
		double err_x=0;
		double err_y=0;
		double err_xx=0;
		double err_xy=0;
		double err_yy=0;

		// loop over xi_1 and xi_2, i.e. left/right sides of cartesian grid
		for (int j=0;j<gr_u.n;j++) {

			// element index
			int i = left ? 0 : (gr_u.m-1);

			// left/right elements
			element eu = gr_u(i,j);

			// trafos
			ct_u = dm->alloc_element_geometry(eu);

			// grab vel solution
			for (node nu:eu) {
				usol[nu.l] = ux[nu.g];
				vsol[nu.l] = uy[nu.g];
			}

			// element coords
			double xi = (left?0:1),et;

			// loop through gauss points
			for (int g=0;g<q.n;g++) {

				// we're integrating over the eta value
				et = q.x[g];
				
				// get Jacobean info
				ct_u->Jinv(xi,et,xi_x,et_x,xi_y,et_y);
				ct_u->xi2x(xi,et,x,y);
				double tmp_det = xi_x*et_y - et_x*xi_y;

				// get polar value
				double th = atan2(y,x-(left?xl:xr));

				// construct line segment length
				x_et =-xi_y/tmp_det;
				y_et = xi_x/tmp_det; 

				// relevant determinant is 1D in xi-direction
				det = sqrt(x_et*x_et + y_et*y_et) * fabs(y) * M_PI;
				
				// construct normal unit vector
				// from xi-gradient
	//			nx = xi_x; ny = xi_y;

				// grab normal from position
				double nnx = x - (left?xl:xr);
				double nny = y;
				double nnmag = sqrt(nnx*nnx + nny*nny);
				nnx /= nnmag; nny /= nnmag;
				nx = nnx; ny = nny;


	//			if (!left) {nx *= -1; ny *= -1;}
				
				// normalize
	//			n2 = nx*nx + ny*ny;
	//			nx /= sqrt(n2);
	//			ny /= sqrt(n2);

	//			mesg("%g (%g), %g (%g)\n",nx,nnx,ny,nny);

//				nx = nnx; ny = nny;

				// normal derivative
				// (only this special form here!)
//				nx_x =  ny*ny / R;
//				ny_x = -nx*ny / R;

				// pull out vel and pressure sol'ns
				mast_u.eval<dcomplex>(*ct_u,xi,et,usol,ux_v,ux_x,ux_y,ux_xx,ux_xy,ux_yy);
				mast_u.eval<dcomplex>(*ct_u,xi,et,vsol,uy_v,uy_x,uy_y,uy_xx,uy_xy,uy_yy);

				brink_derivs(x,y,
						ux_k_v,ux_k_x,ux_k_y,ux_k_xx,ux_k_xy,ux_k_yy,
						uy_k_v,uy_k_x,uy_k_y,uy_k_xx,uy_k_xy,uy_k_yy);

//				mesg("x = %10.6g, y = %10.6g, ux = %10.6g + i * %10.6g, ux_k = %10.6g + i * %10.6g\n",
//					x,y,std::real(ux_v),std::imag(ux_v),std::real(ux_k_v),std::imag(ux_k_v));

				ux_v -= ux_k_v;
				ux_x -= ux_k_x;
				ux_y -= ux_k_y;
				ux_xx -= ux_k_xx;
				ux_xy -= ux_k_xy;
				ux_yy -= ux_k_yy;

				uy_v -= uy_k_v;
				uy_x -= uy_k_x;
				uy_y -= uy_k_y;
				uy_xx -= uy_k_xx;
				uy_xy -= uy_k_xy;
				uy_yy -= uy_k_yy;

				err_v += std::real(ux_v * std::conj(ux_v) + uy_v * std::conj(uy_v)) * q.w[g] * det;
				err_x += std::real(ux_x * std::conj(ux_x) + uy_x * std::conj(uy_x)) * q.w[g] * det;
				err_y += std::real(ux_y * std::conj(ux_y) + uy_y * std::conj(uy_y)) * q.w[g] * det;
				err_xx += std::real(ux_xx * std::conj(ux_xx) + uy_xx * std::conj(uy_xx)) * q.w[g] * det;
				err_xy += std::real(ux_xy * std::conj(ux_xy) + uy_xy * std::conj(uy_xy)) * q.w[g] * det;
				err_yy += std::real(ux_yy * std::conj(ux_yy) + uy_yy * std::conj(uy_yy)) * q.w[g] * det;



			}

			dm->free_element_geometry(ct_u);
		}

		delete[] vsol;
		delete[] usol;

		if (gru()) {
			fprintf(fh,"%d %g %g %g %g %g %g\n",m,
				sqrt(err_v),
				sqrt(err_x),
				sqrt(err_y),
				sqrt(err_xx),
				sqrt(err_xy),
				sqrt(err_yy));
		}

		/*
		mesg("(ux,        uxx,      p,     px): %10.6g %10.6g %10.6g %10.6g\n",std::real(Fux),std::real(Fuxx),std::real(Fp),std::real(Fpx));
		mesg("(uy,        uxy,     vx,    vxx): %10.6g %10.6g %10.6g %10.6g\n",std::real(Fuy),std::real(Fuxy),std::real(Fvx),std::real(Fvxx));
		mesg("(ds.n x, ds.n y, s.dn x, s.dn y): %10.6g %10.6g %10.6g %10.6g\n",std::real(F1),std::real(F2),std::real(F3),std::real(F4));
		// */
		
	}

	void calc_error(FILE *fh) {

//		mesg("boop\n");

		// arrays for mass matrix
		double *mm_u = new double[lu*lu];
		double *mm_p = new double[lp*lp];

		// local arrays for solution
		dcomplex *usol = new dcomplex[lu],
			     *vsol = new dcomplex[lu],
				 *psol = new dcomplex[lp];

		dcomplex *du = new dcomplex[lu],
				 *dv = new dcomplex[lu];

		dcomplex *du_x = new dcomplex[lu],
				 *dv_x = new dcomplex[lu];

		dcomplex *du_y = new dcomplex[lu],
				 *dv_y = new dcomplex[lu];

		dcomplex *du_xx = new dcomplex[lu],
				 *dv_xx = new dcomplex[lu];
		dcomplex *du_xy = new dcomplex[lu],
				 *dv_xy = new dcomplex[lu];
		dcomplex *du_yy = new dcomplex[lu],
				 *dv_yy = new dcomplex[lu];

		// local arrays for errors
		double *verr = new double[lu],
		       *perr = new double[lp];

		// variables for solutions
		dcomplex uu,uu_x,uu_y,uu_xx,uu_xy,uu_yy;
		dcomplex vv,vv_x,vv_y,vv_xx,vv_xy,vv_yy;
		dcomplex pp,pp_x,pp_y,pp_xx,pp_xy,pp_yy;

		trafo_base *ct_u,*ct_p;

		double verr_v_tot = 0;
		double verr_x_tot = 0;
		double verr_y_tot = 0;
		double verr_xx_tot = 0;
		double verr_xy_tot = 0;
		double verr_yy_tot = 0;

		for (element eu:gr_u) {

//			if (eu.i == 0 || eu.i == m-1) continue;

			element ep = gr_p(eu.i,eu.j);

			// grab geometry
			ct_u = dm->alloc_element_geometry(eu);
			ct_p = dm->alloc_element_geometry(ep);

			// loop over quadrature points
			for (int qi=0; qi<mast_u.g->N; qi++) {

				// coords
				double xi,et,x,y;
				mast_u.g->point(qi,xi,et);
				ct_u->xi2x(xi,et,x,y);

				// weighting + determinant
				double pre = mast_u.g->weight(qi) * ct_u->Jdet(xi,et);

				// evaluate velocities
				dcomplex ux_v, ux_x, ux_y, ux_xx, ux_xy, ux_yy;
				dcomplex uy_v, uy_x, uy_y, uy_xx, uy_xy, uy_yy;

				// evaluate known velocities
				dcomplex ux_k_v, ux_k_x, ux_k_y, ux_k_xx, ux_k_xy, ux_k_yy;
				dcomplex uy_k_v, uy_k_x, uy_k_y, uy_k_xx, uy_k_xy, uy_k_yy;

				// grab solutions 
				for (node n:eu) {

					// get vel value at node
					usol[n.l] = ux[n.g];

//					printf("li=%d, lg=%d, v=%g + i*%g (Vl = %g + i*%g, Vr = %g + i*%g)\n",n.l,n.g,ux[n.g],
//							std::real(V_l.x),std::imag(V_l.x),std::real(V_r.x),std::imag(V_r.x));

					vsol[n.l] = uy[n.g];
				}

				mast_u.eval(*ct_u,xi,et,usol,ux_v,ux_x,ux_y,ux_xx,ux_xy,ux_yy,true);
				mast_u.eval(*ct_u,xi,et,vsol,uy_v,uy_x,uy_y,uy_xx,uy_xy,uy_yy);

				brink_derivs(x,y,ux_k_v,ux_k_x,ux_k_y,ux_k_xx,ux_k_xy,ux_k_yy,
					uy_k_v,uy_k_x,uy_k_y,uy_k_xx,uy_k_xy,uy_k_yy);

				if (std::isnan(std::real(ux_k_v)) || std::isnan(std::imag(ux_k_v))) mesg("nan velocity at x=%g,y=%g (xl=%g,xr=%g)\n",x,y,xl,xr);
				if (std::isnan(std::real(ux_k_x)) || std::isnan(std::imag(ux_k_x))) mesg("nan velocity at x=%g,y=%g (xl=%g,xr=%g)\n",x,y,xl,xr);

//				printf("ux=%g, uy=%g, ux*=%g, uy*=%g\n",
//						std::real(ux_v), std::real(uy_v),
//						std::real(ux_k_v), std::real(uy_k_v));

				
//				mesg("x = %10.6g, y = %10.6g, ux_xx = %10.6g + i * %10.6g, ux_k_xx = %10.6g + i * %10.6g\n",
//					x,y,std::real(ux_xx),std::imag(ux_xx),std::real(ux_k_xx),std::imag(ux_k_xx));


				// subtract off
				ux_v -= ux_k_v;
				ux_x -= ux_k_x;
				ux_y -= ux_k_y;
				ux_xx -= ux_k_xx;
				ux_xy -= ux_k_xy;
				ux_yy -= ux_k_yy;

				uy_v -= uy_k_v;
				uy_x -= uy_k_x;
				uy_y -= uy_k_y;
				uy_xx -= uy_k_xx;
				uy_xy -= uy_k_xy;
				uy_yy -= uy_k_yy;


				// add to totals

				verr_v_tot += std::real((ux_v*std::conj(ux_v) + uy_v*std::conj(uy_v))) * pre;
				verr_x_tot += std::real((ux_x*std::conj(ux_x) + uy_x*std::conj(uy_x))) * pre;
				verr_y_tot += std::real((ux_y*std::conj(ux_y) + uy_y*std::conj(uy_y))) * pre;
				verr_xx_tot += std::real((ux_xx*std::conj(ux_xx) + uy_xx*std::conj(uy_xx))) * pre;
				verr_xy_tot += std::real((ux_xy*std::conj(ux_xy) + uy_xy*std::conj(uy_xy))) * pre;
				verr_yy_tot += std::real((ux_yy*std::conj(ux_yy) + uy_yy*std::conj(uy_yy))) * pre;
			}

			/*
			// grab mass matrices
			mast_u.mass_matrix(*ct_u,mm_u);
			mast_p.mass_matrix(*ct_p,mm_p);

			// get velocity solution
			for (node n:eu) {

				// get vel value at node
				usol[n.l] = ux[n.g];
				vsol[n.l] = uy[n.g];

				// get cartesian coord
				double xi,et,x,y;
				mast_u.point(n.l,xi,et);
				ct_u->xi2x(xi,et,x,y);

				// get exact solution
				vcomplex vc = brink_sol(x,y);

				// get pieces of error
				du[n.l] = ux[n.g] - vc.x;
				dv[n.l] = uy[n.g] - vc.y;

//				mesg("element (%d,%d), node (%d,%d): x=%g, y=%g, u = %g + %g*j, v = %g + %g*j, u* = %g + %g*j, v* = %g + %g*j\n",
//						x,y,
//					eu.i,eu.j,n.li,n.lj, std::real(ux[n.g]), std::imag(ux[n.g]), std::real(uy[n.g]), std::imag(uy[n.g]),
//					std::real(vc.x), std::imag(vc.x), std::real(vc.y), std::imag(vc.y));

//				mesg("element (%d,%d), node (%d,%d): du = %g + %g*j, dv = %g + %g*j\n",
//					eu.i,eu.j,n.li,n.lj, std::real(du[n.l]), std::imag(du[n.l]), std::real(dv[n.l]), std::imag(dv[n.l]));

//				mesg("%d: %g\n",n.l,mm_u[n.l + lu*n.l]);
			}

			// get velocity derivatives 
			for (node n:eu) {

				// get cartesian coord
				double xi,et,x,y;
				mast_u.point(n.l,xi,et);

				dcomplex a_,b_,c_,d_,e_;
				mast_u.eval<dcomplex>(*ct_u,xi,et,usol,a_, du_x[n.l],
						du_y[n.l],du_xx[n.l],du_xy[n.l],du_yy[n.l]);
				mast_u.eval<dcomplex>(*ct_u,xi,et,vsol,a_, dv_x[n.l],
						dv_y[n.l],dv_xx[n.l],dv_xy[n.l],dv_yy[n.l]);

				ct_u->xi2x(xi,et,x,y);

				// get exact solution
				dcomplex ux_x,ux_y,ux_xx,ux_xy,ux_yy;
				dcomplex uy_x,uy_y,uy_xx,uy_xy,uy_yy;
				brink_derivs(x,y,
					ux_x,ux_y,ux_xx,ux_xy,ux_yy,
					uy_x,uy_y,uy_xx,uy_xy,uy_yy);

				// get pieces of error
				du_x[n.l] -= ux_x;
				dv_x[n.l] -= uy_x;
				du_y[n.l] -= ux_y;
				dv_y[n.l] -= uy_y;
				du_xx[n.l] -= ux_xx;
				dv_xx[n.l] -= uy_xx;
				du_xy[n.l] -= ux_xy;
				dv_xy[n.l] -= uy_xy;
				du_yy[n.l] -= ux_yy;
				dv_yy[n.l] -= uy_yy;

			}

			for (node nj:eu) for (node ni:eu) {

				dcomplex ui(du[ni.l]), uj(du[nj.l]);
				dcomplex vi(dv[ni.l]), vj(dv[nj.l]);

				dcomplex ui_x(du_x[ni.l]), uj_x(du_x[nj.l]);
				dcomplex vi_x(dv_x[ni.l]), vj_x(dv_x[nj.l]);

				dcomplex ui_y(du_y[ni.l]), uj_y(du_y[nj.l]);
				dcomplex vi_y(dv_y[ni.l]), vj_y(dv_y[nj.l]);

				dcomplex ui_xx(du_xx[ni.l]), uj_xx(du_xx[nj.l]);
				dcomplex vi_xx(dv_xx[ni.l]), vj_xx(dv_xx[nj.l]);

				dcomplex ui_xy(du_xy[ni.l]), uj_xy(du_xy[nj.l]);
				dcomplex vi_xy(dv_xy[ni.l]), vj_xy(dv_xy[nj.l]);

				dcomplex ui_yy(du_yy[ni.l]), uj_yy(du_yy[nj.l]);
				dcomplex vi_yy(dv_yy[ni.l]), vj_yy(dv_yy[nj.l]);

//				dcomplex ui(1-I),uj(1-I);
//				dcomplex vi(1-I),vj(1-I);

				verr_tot +=
					(ui * std::conj(uj) + vi * std::conj(vj)) * mm_u[ni.l + lu * nj.l];

				verr_x_tot +=
					(ui_x * std::conj(uj_x) + vi_x * std::conj(vj_x)) * mm_u[ni.l + lu * nj.l];
				verr_y_tot +=
					(ui_y * std::conj(uj_y) + vi_y * std::conj(vj_y)) * mm_u[ni.l + lu * nj.l];
				verr_xx_tot +=
					(ui_xx * std::conj(uj_xx) + vi_xx * std::conj(vj_xx)) * mm_u[ni.l + lu * nj.l];
				verr_xy_tot +=
					(ui_xy * std::conj(uj_xy) + vi_xy * std::conj(vj_xy)) * mm_u[ni.l + lu * nj.l];
				verr_yy_tot +=
					(ui_yy * std::conj(uj_yy) + vi_yy * std::conj(vj_yy)) * mm_u[ni.l + lu * nj.l];

//				mesg("%d %d: %g %g\n",ni.l,nj.l,mm_u[ni.l + lu*nj.l],mm_u[nj.l + lu*ni.l]);

//				mesg("%g %g\n",std::real(verr_tot),std::imag(verr_tot));

			}
			*/

			dm->free_element_geometry(ct_u);
			dm->free_element_geometry(ct_p);
//			mesg("%g %g\n",std::real(verr_tot),std::imag(verr_tot));
		}


		if (gru()) {

			fprintf(fh,"%d %g %g %g %g %g %g\n",m,
			sqrt(verr_v_tot),
			sqrt(verr_x_tot),
			sqrt(verr_y_tot),
			sqrt(verr_xx_tot),
			sqrt(verr_xy_tot),
			sqrt(verr_yy_tot));

			fflush(fh);
		}

		// cleanup
		delete[] du;
		delete[] dv;

		delete[] perr;
		delete[] verr;

		delete[] usol;
		delete[] vsol;
		delete[] psol;

		delete[] mm_p;
		delete[] mm_u;

//		mesg("boop bop\n");
	}

	// secant method force condition solves
	void solve_force_balance(dcomplex Fmag=1);
	void solve_force_free(bool approx=false);
	void solve_force_osc();

	// coarsening
	brinkman *coarsen(int cx,int cy);
};

// shorthand for pair of C-strings
typedef std::pair<const char*,const char*> string_duple;


#endif
