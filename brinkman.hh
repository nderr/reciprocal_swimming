///////////////////////////////////////////////////////////////////////////////
//
// brinkman.hh
//
// a class representing a complex-valued Brinkman equation:
// 
// (nabla^2 - alpha^2) u = grad p - f
//
///////////////////////////////////////////////////////////////////////////////

#ifndef BRINKMAN_HH
#define BRINKMAN_HH

#include <cassert>
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
	/** shortcut for zero vector */
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
	// calculate the total problem matrix
	static void calculate_A_sys(DM da_sys,Mat P,brinkman*);

	// caclulate the RHS vector
	static PetscErrorCode calculate_B(KSP ksp,Vec B,void *ctx);

	/**
	 * Verifies PETSc solve completed successfully
	 */
	void check_sol() {

		KSPConvergedReason reason;
		KSPGetConvergedReason(ksp,&reason);

		if (reason>0) return;
		
		KSPDestroy(&ksp);
		bail(-1,"didn't converge for reason %d\n",reason);
	}

	/**
	 * helper function to sum a vector of integers
	 */
	static int sum(std::vector<int> v) {
		return std::accumulate(v.begin(),v.end(),0);
	}

	/**
	 * helper function to calculate numbers of processors
	 */
	static void proc_counts(brinkman *prob,
		std::vector<int> &lux,std::vector<int> &luy,
		std::vector<int> &lpx,std::vector<int> &lpy);

	/**
	 * helper function to find best processor decomposition
	 */
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

	/** helper for setting <Vx-Vr> = 1 */
	void norm_Vx() {

		// get magnitude of oscillation magnitude
		dcomplex sum(V_r.x-V_l.x);

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

	dcomplex avg_vel(bool left);
	/* force-free integration functions */
	dcomplex drag_force(bool left,bool approx=false);
	/* calculate the exact or approximate drag coefficient */
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

	// secant method force condition solves
	void solve_force_balance(dcomplex Fmag=1);
	void solve_force_free();

	// coarsening
	brinkman *coarsen(int cx,int cy);
};

// shorthand for pair of C-strings
typedef std::pair<const char*,const char*> string_duple;

#endif
