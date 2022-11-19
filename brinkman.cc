//////////////////////////////////////////////////////////////////////////////
//
// brinkman.c
//
// implementations for solving the Brinkman equation
// 
// (nabla^2 - alpha^2) u = grad p - f
//
///////////////////////////////////////////////////////////////////////////////

#include "brinkman.hh"
#include <unistd.h>

//////////////////// top-level and body-force ////////////////////

/**
 * helper function returning whether we're on zero proc
 *
 * @return if we're on zero proc
 */
bool gru() {
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	return rank==0;
}

/**
 * helper function for bailing
 *
 * @param[in] code exit code
 * @param[in] fmt printf-style format code
 * @param[in] ... formatted args
 */
void bail(int code,const char *fmt,...) {

	if (gru()) {
		// set up variable args
		va_list args; va_start(args,fmt);
		// print and finish
		vfprintf(stderr,fmt,args); va_end(args);
	}

	MPI_Barrier(PETSC_COMM_WORLD);
	PetscFinalize();

	// exit
	exit(code);
}

/**
 * helper function for printing on master proc
 *
 * @param[in] fmt printf-style format code
 * @param[in] ... formatted args
 */
void mesg(const char *fmt,...) {

	// only print on proc 0
	if (!gru()) return;

	// set up variable args
	va_list args; va_start(args,fmt);
	// print and finish
	vfprintf(stdout,fmt,args); va_end(args);
}

///////////////////// vcomplex /////////////////////////

/** complex zero vector */
vcomplex vcomplex::zero(0);


/**
 * body-force source function at a given location representing
 * a Reynolds stress due to a leading-order Brinkman solution
 *
 * @param[in] x x-coordinate
 * @param[in] y y-coordinate
 * @return body-force at location (x,y)
 */
vcomplex reynolds::operator()(double x,double y) const {

	double Re = br->alpha * br->alpha;

	// body force (Reynolds stress)
	vcomplex f = 0.5 * Re * br->u_grad_u(x,y);

	// return real part
	return f.real();
}

/////////////////// brinkman //////////////////////////

/**
 * General constructor for brinkman
 */
brinkman::brinkman(ptype pt,int po_,int m_,int n_,
	bool hp_,bool vp_,domain *dm_,tsrc_base<vcomplex> *f_,
	double mu_,double a_,vcomplex V):
	type(pt),po(po_),m(m_),n(n_),hp(hp_),vp(vp_),
	gr_u(po,m,n,hp,vp),gr_p(po-1,m,n,hp,vp),
	mast_u(gr_u.p),    mast_p(gr_p.p,gr_u.p),
	ump(gr_u.mp),unp(gr_u.np),lu((1+po)*(1+po)),gu(ump*unp),
	pmp(gr_p.mp),pnp(gr_p.np),lp(po*po),        gp(pmp*pnp),
	dm(dm_),rev(dynamic_cast<cylindrical*>(dm)==NULL?false:true),
	f(f_),mu(mu_),alpha(a_),V_inf(V),data(false),ux(NULL),uy(NULL),p(NULL){}

/**
 * General destructor
 */
brinkman::~brinkman() {
	if (data) {
		delete[] ux;
		delete[] uy;
		delete[] p;
	}
}

/**
 * allocates space for the solution to a Brinkman equation
 */
void brinkman::alloc_data() {
	if (!data) {
		data=true;
		ux = new dcomplex[gu];
		uy = new dcomplex[gu];
		p  = new dcomplex[gp];
	}
}

/**
 * writes data to an open binary file
 *
 * @param[in] fb binary FILE to write to
 */
void brinkman::write_data(FILE *fb) {
	fwrite(ux,sizeof(dcomplex),gu,fb);
	fwrite(uy,sizeof(dcomplex),gu,fb);
	fwrite(p ,sizeof(dcomplex),gp,fb);
}

/**
 * reads data from an open binary file
 *
 * @param[in] fb binary FILE to read from
 */
void brinkman::read_data(FILE *fb) {
	alloc_data();
	fread(ux,sizeof(dcomplex),gu,fb);
	fread(uy,sizeof(dcomplex),gu,fb);
	fread(p ,sizeof(dcomplex),gp,fb);
}

/**
 * Static function for problem construction. Takes in a DMDA
 * describing a nodal grid of specified polynomial order and
 * returns the range of ownership of finite elements. It's
 * critical for matrix / source vector computation that each finite
 * element in the grid be represented on only one processor
 *
 * @param[in] dmda array of nodes
 * @param[in] poly polynomial order
 * @param[out] eia first (inclusive) x-element on this proc
 * @param[out] eib first (exclusive) x-element NOT on this proc
 * @param[out] eja first (inclusive) y-element on this proc
 * @param[out] ejb first (exclusive) y-element NOT on this proc
 */
void brinkman::dm_fe_range(DM dmda,int p,
		int &eia,int &eib,int &eja,int &ejb) {

	// grab grid info
	DMDALocalInfo info; DMDAGetLocalInfo(dmda,&info);

	// what is grid size (node-wise)
	int mp(info.mx),np(info.my);

	// if not periodic, subtract off endpoint
	if (info.bx!=DM_BOUNDARY_PERIODIC) mp--;
	if (info.by!=DM_BOUNDARY_PERIODIC) np--;

	// divide by poly order for num elements
	int m(mp/p),n(np/p);

	eia=info.xs/p; eib=(info.xs+info.xm-1)/p + 1; if(eib>m)eib=m;
	eja=info.ys/p; ejb=(info.ys+info.ym-1)/p + 1; if(ejb>n)ejb=n;
}


/**
 * Static function which takes in a node, element, DMDALocalInfo struct,
 * and component index if order to calculate the local processor indices
 * of the node's provided component
 *
 * @param[in] n node in grid
 * @param[in] e finite element where node n is based
 * @param[in] info DMDALocalInfo corresponding to grid e and n are in
 * @param[in] c component index
 */
int brinkman::node_local_ind(node &n,element &e,DMDALocalInfo &info,int c) {

	// check component represents valid field
	assert(c >= 0 && c < info.dof);

	// grab global index (grab ghost row if at top of periodic grid)
	int gi = (e.h_prd && n.gi==0 && n.li!=0) ? (e.p*e.m) : n.gi;
	int gj = (e.v_prd && n.gj==0 && n.lj!=0) ? (e.p*e.n) : n.gj;

	// ordering is components, then x, then y, etc
	return ((gi-info.gxs) + (gj-info.gys)*info.gxm)*info.dof + c;
}

/**
 * Static function which takes in a node, element, DMDALocalInfo struct,
 * and component index if order to calculate the local processor indices
 * of the node's provided component
 *
 * @param[in] n node in grid
 * @param[in] e finite element where node n is based
 * @param[in] info DMDALocalInfo corresponding to grid e and n are in
 * @param[in] c component index
 */
int brinkman::node_natural_ind(node &n,element &e,DMDALocalInfo &info,int c) {

	// check component represents valid field
	assert(c >= 0 && c < info.dof);

	// stride len
	int ml = e.h_prd ? (e.m*e.p) : (e.m*e.p+1);

	// ordering is components, then x, then y, etc
	return (n.gi + n.gj*ml)*info.dof + c;
}

/**
 * Get local indices (i.e. indices corresponding to position on current
 * processor) for the x-velocity field in the provided element
 *
 * @param[in] e element
 * @param[in] da_vel velocity nodes DMDA
 * @param[out] ind the velocity indices
 */
void brinkman::get_u_local_indices(element &e,DM da_vel,int *ind) {

	// grab grid info
	DMDALocalInfo info;
	DMDAGetLocalInfo(da_vel,&info);

	// fill in with local indices
	for (node n:e) ind[n.l] = node_local_ind(n,e,info,0);
}

/**
 * Get local indices (i.e. indices corresponding to position on current
 * processor) for the y-velocity field in the provided element
 *
 * @param[in] e element
 * @param[in] da_vel velocity nodes DMDA
 * @param[out] ind the velocity indices
 */
void brinkman::get_v_local_indices(element &e,DM da_vel,int *ind) {
	DMDALocalInfo info;
	DMDAGetLocalInfo(da_vel,&info);
	for (node n:e) ind[n.l] = node_local_ind(n,e,info,1);
}

/**
 * Get local indices (i.e. indices corresponding to position on current
 * processor) for the pressure field in the provided element
 *
 * @param[in] e element
 * @param[in] da_pre pressure nodes DMDA
 * @param[out] ind the pressure indices
 */
void brinkman::get_p_local_indices(element &e,DM da_pre,int *ind) {
	DMDALocalInfo info;
	DMDAGetLocalInfo(da_pre,&info);
	for (node n:e) ind[n.l] = node_local_ind(n,e,info);
}

/**
 * Get natural indices (i.e. indices corresponding to position w/in
 * global node grid) for the x-velocity field in the provided element
 *
 * @param[in] e element
 * @param[in] da_vel velocity nodes DMDA
 * @param[out] ind the velocity indices
 */
void brinkman::get_u_natural_indices(element &e,DM da_vel,int *ind) {

	// grab grid info
	DMDALocalInfo info;
	DMDAGetLocalInfo(da_vel,&info);

	// fill in with natural indices
	for (node n:e) ind[n.l] = node_natural_ind(n,e,info,0);
}

/**
 * Get natural indices (i.e. indices corresponding to position w/in
 * global node grid) for the y-velocity field in the provided element
 *
 * @param[in] e element
 * @param[in] da_vel velocity nodes DMDA
 * @param[out] ind the velocity indices
 */
void brinkman::get_v_natural_indices(element &e,DM da_vel,int *ind) {
	DMDALocalInfo info;
	DMDAGetLocalInfo(da_vel,&info);
	for (node n:e) ind[n.l] = node_natural_ind(n,e,info,1);
}

/**
 * Get natural indices (i.e. indices corresponding to position w/in
 * global node grid) for the pressure field in the provided element
 *
 * @param[in] e element
 * @param[in] da_pre pressure nodes DMDA
 * @param[out] ind the pressure indices
 */
void brinkman::get_p_natural_indices(element &e,DM da_pre,int *ind) {
	DMDALocalInfo info;
	DMDAGetLocalInfo(da_pre,&info);
	for (node n:e) ind[n.l] = node_natural_ind(n,e,info);
}

/**
 * Get global indices (i.e. indices corresponding to position w/in
 * concatenation of each proc's local node ordering) for the x-velocity
 * field in the provided element
 *
 * @param[in] e element
 * @param[in] da_vel velocity nodes DMDA
 * @param[out] ind the velocity indices
 */
void brinkman::get_u_global_indices(element &e,DM da_vel,int *ind) {
	get_u_natural_indices(e,da_vel,ind);
	AO ao;
	DMDAGetAO(da_vel,&ao);
	AOApplicationToPetsc(ao,lu,ind);
}

/**
 * Get global indices (i.e. indices corresponding to position w/in
 * concatenation of each proc's local node ordering) for the y-velocity
 * field in the provided element
 *
 * @param[in] e element
 * @param[in] da_vel velocity nodes DMDA
 * @param[out] ind the velocity indices
 */
void brinkman::get_v_global_indices(element &e,DM da_vel,int *ind) {
	get_v_natural_indices(e,da_vel,ind);
	AO ao;
	DMDAGetAO(da_vel,&ao);
	AOApplicationToPetsc(ao,lu,ind);
}

/**
 * Get global indices (i.e. indices corresponding to position w/in
 * concatenation of each proc's local node ordering) for the pressure
 * field in the provided element
 *
 * @param[in] e element
 * @param[in] da_pre pressure nodes DMDA
 * @param[out] ind the pressure indices
 */
void brinkman::get_p_global_indices(element &e,DM da_pre,int *ind) {
	get_p_natural_indices(e,da_pre,ind);
	AO ao;
	DMDAGetAO(da_pre,&ao);
	AOApplicationToPetsc(ao,lp,ind);
}

/**
 * Calculate the source vector components corresponding to the
 * integral over the provided finite element
 *
 * @param[in] eu finite element
 * @param[out] fx x-directed body force source
 * @param[out] fy y-directed body force source
 */
void brinkman::get_BU(element &eu,PetscScalar *fx,PetscScalar *fy) {
	
	// grab element geometry and source componenets
	trafo_base *ct = dm->alloc_element_geometry(eu);
	std::vector<vcomplex> vc(lu);
	mast_u.src_term(*ct,*f,vc.data());

	// coppy into output arrays
	PetscScalar *xp(fx),*yp(fy);
	for(vcomplex v:vc) { (*fx++)=v.x; (*fy++)=v.y; }

	// clean up
	dm->free_element_geometry(ct);
}

/**
 * get u-u part of problem matrix for provided element
 *
 * @param[in] eu element in which to get matrix
 * @param[out] m local u-u part of matrix
 */
void brinkman::get_Auu(element &eu,PetscScalar *m) {

	// grab element geometry, double mat
	trafo_base *ct = dm->alloc_element_geometry(eu);
	double *v = new double[lu*lu];

	// get stiffness matrix, copy into scalars
	mast_u.stiffness_matrix(*ct,v);

	// weight by viscosity
	for(int i=0;i<lu*lu;i++) m[i] = mu * v[i];

	// get mass matrix
	mast_u.mass_matrix(*ct,v);

	// weight by imaginary Reynolds
	const double Re = alpha * alpha;
	for(int i=0;i<lu*lu;i++) m[i] += I*Re*v[i];

	// clean up
	delete[] v;
	dm->free_element_geometry(ct);
}

/**
 * get v-v part of problem matrix for provided element
 *
 * @param[in] eu element in which to get matrix
 * @param[out] m local v-v part of matrix
 */
void brinkman::get_Avv(element &eu,PetscScalar *m) {

	// start with x-mat
	get_Auu(eu,m);

	// only do more if revolved
	if (!rev) return;

	// grab element geometry, double mat
	trafo_base *ct = dm->alloc_element_geometry(eu);
	double *v = new double[lu*lu];

	// get revolved mass mat, copy into scalars
	mast_u.mass_matrix_cyl_y(*ct,v); // third arg -> add vs insert

	// weight by viscosity
	for(int i=0;i<lu*lu;i++) m[i] += mu * v[i];

	// clean up
	delete[] v;
	dm->free_element_geometry(ct);
}

/**
 * get u-p and v-p part of problem matrix for provided element
 *
 * @param[in] eu element in which to get matrix
 * @param[out] mx local u-p part of matrix
 * @param[out] my local v-p part of matrix
 * @param[in] T whether to get transposes
 */
void brinkman::get_Auvp(element &eu,PetscScalar *mx,PetscScalar *my,bool T) {

	// grab element geometry
	trafo_base *ct = dm->alloc_element_geometry(eu);
	
	// double mats for master values
	double *vx = new double[lu*lp], *vy = new double[lu*lp];

	// get b-matrices
	master::stokes_bmats(mast_u,mast_p,*ct,vx,vy,rev);

	// copy in -> set BOTH NEGATIVE
	for (int r=0,in=0;r<lu;r++) {
		for (int c=0;c<lp;c++,in++) {

			int out = T ? (r+c*lu) : (c+r*lp);

			mx[out] = -vx[in];
			my[out] = -vy[in];
		}
	}

	// clean up
	delete[] vy;
	delete[] vx;
	dm->free_element_geometry(ct);
}

/**
 * get this element's part of the pressure mass matrix
 *
 * @param[in] ep element in which to get matrix
 * @param[out] m local p-p part of matrix
 */
void brinkman::get_P_mass(element &ep,PetscScalar *m) {

	// grab element geometry, double mat
	trafo_base *ct = dm->alloc_element_geometry(ep);
	double *v = new double[lp*lp];

	// get mass matrix
	mast_p.mass_matrix(*ct,v);

	// weight by inverse viscosity
	for(int i=0;i<lp*lp;i++) m[i] = v[i] / mu;

	// clean up
	delete[] v;
	dm->free_element_geometry(ct);
}

/**
 * get this element's part of the pressure stiffness matrix
 *
 * @param[in] ep element in which to get matrix
 * @param[out] m local p-p part of matrix
 */
void brinkman::get_P_stiff(element &ep,PetscScalar *m) {

	// grab element geometry, double mat
	trafo_base *ct = dm->alloc_element_geometry(ep);
	double *v = new double[lp*lp];

	// get stiffness matrix
	mast_p.stiffness_matrix(*ct,v);

	// weight by inverse drag length
	for(int i=0;i<lp*lp;i++) m[i] = v[i] / (alpha*alpha*I);

	// clean up
	delete[] v;
	dm->free_element_geometry(ct);
}

/**
 * template function for retrieving and inserting element-wise
 * portions of velocity problem matrix, with rows and columns
 * left unassembled at nodes representing Dirichlet BCs
 *
 * @tparam x_vel (otherwise y vel)
 * @param[in] eu element in which to get matrix
 * @param[in] da_vel velocity node DMDA
 * @param[out] Auu local velocity-velocity submatrix
 */
template<bool x_vel>
void brinkman::add_AUU(element &eu,DM da_vel,Mat Auu) {

	// non-fixed indices
	std::vector<int> iu_nf;

	// non-fixed mat and source values
	std::vector<dcomplex> m_nf;

	// get matrix values
	std::vector<dcomplex> m_all(lu*lu);
	if (x_vel) get_Auu(eu,m_all.data());
	else       get_Avv(eu,m_all.data());

	// get local indices
	std::vector<int> li(lu);
	if (x_vel) get_u_local_indices(eu,da_vel,li.data());
	else       get_v_local_indices(eu,da_vel,li.data());

	// pointers to functions about fixed points
	bool (brinkman::*fix_pt)(node &n) = 
		x_vel ? &brinkman::fix_pt_ux : &brinkman::fix_pt_uy;
	dcomplex (brinkman::*fix_vl)(node &n) =
		x_vel ? &brinkman::fix_vl_ux : &brinkman::fix_vl_uy;

	// loop through non-fixed (row) nodes
	// store index/vals as appropriate
	for(node nr:eu) if (!(this->*fix_pt)(nr)) {

		// store non-fixed index
		iu_nf.push_back(li[nr.l]);

		// store each non-fixed matrix column
		for(node nc:eu) if (!(this->*fix_pt)(nc))
			 m_nf.push_back(m_all[nc.l+nr.l*lu]);
	}

	// add non-fixed values
	int ni(iu_nf.size()),*ii(iu_nf.data());
	dcomplex *d(m_nf.data());
	MatSetValuesLocal(Auu,ni,ii,ni,ii,d,ADD_VALUES);
}

/**
 * function for retrieving and inserting element-wise
 * portions of the off-diagonal velocity-pressure problem
 * matrix, with rows and columns left unassembled at nodes
 * representing Dirichlet BCs
 *
 * @param[in] eu velocity element in which to get matrix
 * @param[in] ep pressure element in which to get matrix
 * @param[in] da_vel velocity node DMDA
 * @param[in] da_pre pressure node DMDA
 * @param[out] Aup local velocity-pressure submatrix
 * @param[out] Apu local pressure-velocity submatrix
 */
void brinkman::add_Auvp(element &eu,element &ep,
		DM da_vel,DM da_pre,Mat Aup,Mat Apu) {

	// non-fixed vel indices
	std::vector<int> i_nf;

	// non-fixed mat and source values
	std::vector<dcomplex> m_nf;

	// get matrix values (row vel, col pre)
	dcomplex *mx = new dcomplex[2*lu*lp], *my(mx+lu*lp);
	get_Auvp(eu,mx,my);

	// get local indices
	std::vector<int> liu(lu),liv(lu),lip(lp);
	get_u_local_indices(eu,da_vel,liu.data());
	get_v_local_indices(eu,da_vel,liv.data());
	get_p_local_indices(ep,da_pre,lip.data());

	// loop through vel (row) nodes
	for (node nr:eu) {

		// x-fixing

		// loop through pressure nodes, adding to non-fixed
		// mat or updating source vec as appropriate
		if (!fix_pt_ux(nr)) {
			i_nf.push_back(liu[nr.l]);
			for(node nc:ep) m_nf.push_back(mx[nc.l+nr.l*lp]);
		}

		// y-fixing
		if (!fix_pt_uy(nr)) {
			i_nf.push_back(liv[nr.l]);
			for(node nc:ep) m_nf.push_back(my[nc.l+nr.l*lp]);
		}
	}

	delete[] mx;

	// add non-fixed values
	int np(lip.size()),*ip(lip.data()),nu(i_nf.size()),*iu(i_nf.data());
	dcomplex *d(m_nf.data());
	MatSetValuesLocal(Aup,nu,iu,np,ip,d,ADD_VALUES);

	// transpose the nonzeros
	std::vector<dcomplex> m_nf_T;
	for(size_t c=0;c<lip.size();c++) for(size_t r=0;r<i_nf.size();r++) {
		m_nf_T.push_back(m_nf[c + r*lp]);
	}

	// add non-fixed values
	d = m_nf_T.data();
	MatSetValuesLocal(Apu,np,ip,nu,iu,d,ADD_VALUES);
}

/**
 * function for retrieving and inserting element-wise
 * portions of the pressure mass matrix to be used as
 * a preconditioner for the Schur complement
 *
 * @param[in] e element in which to get matrix
 * @param[in] da_pre pressure node DMDA
 * @param[out] Ppp local pressure-pressure submatrix
 */
void brinkman::add_P_mass(element &ep,DM da_pre,Mat Ppp) {

	// get matrix values
	std::vector<dcomplex> m_all(lp*lp);
	get_P_mass(ep,m_all.data());

	// get local indices
	std::vector<int> li(lp);
	get_p_local_indices(ep,da_pre,li.data());

	// add non-fixed values
	int *ip(li.data());
	dcomplex *d(m_all.data());
	MatSetValuesLocal(Ppp,lp,ip,lp,ip,d,ADD_VALUES);
}

/**
 * function for retrieving and inserting element-wise
 * portions of the pressure mass matrix to be used as
 * a preconditioner for the Schur complement
 *
 * @param[in] e element in which to get matrix
 * @param[in] da_pre pressure node DMDA
 * @param[out] Ppp local pressure-pressure submatrix
 */
void brinkman::add_P_stiff(element &ep,DM da_pre,Mat Ppp) {

	// get matrix values
	std::vector<dcomplex> m_all(lp*lp);
	get_P_stiff(ep,m_all.data());

	// get local indices
	std::vector<int> li(lp);
	get_p_local_indices(ep,da_pre,li.data());

	// add non-fixed values
	int *ip(li.data());
	dcomplex *d(m_all.data());
	MatSetValuesLocal(Ppp,lp,ip,lp,ip,d,ADD_VALUES);
}

/**
 * template function for calculating and inserting portions
 * of velcoity source vector stemming from the zeroing out
 * of columns corresponding to fixed Dirichlet boundary conditions
 *
 * @tparam x_vel (otherwise y vel)
 * @param[in] eu element in which to get matrix
 * @param[in] da_vel velocity node DMDA
 * @param[out] Bu velocity source vector
 */
template<bool x_vel>
void brinkman::add_Diri_BU(element &eu,DM da_vel,Vec Bu) {

	// non-fixed indices
	std::vector<int> iu_nf;

	// non-fixed mat and source values
	std::vector<dcomplex> b;

	// get matrix values
	std::vector<dcomplex> m_all(lu*lu);
	if (x_vel) get_Auu(eu,m_all.data());
	else       get_Avv(eu,m_all.data());

	// get global indices
	std::vector<int> gi(lu);
	if (x_vel) get_u_global_indices(eu,da_vel,gi.data());
	else       get_v_global_indices(eu,da_vel,gi.data());

	// pointers to functions about fixed points
	bool (brinkman::*fix_pt)(node &n) = 
		x_vel ? &brinkman::fix_pt_ux : &brinkman::fix_pt_uy;
	dcomplex (brinkman::*fix_vl)(node &n) =
		x_vel ? &brinkman::fix_vl_ux : &brinkman::fix_vl_uy;

	// loop through non-fixed (row) nodes
	bool adjust = false;
	for(node nr:eu) if (!(this->*fix_pt)(nr)) {

			iu_nf.push_back(gi[nr.l]);
			b.push_back(0);

			// loop through column nodes if unfixed
			for(node nc:eu) if ((this->*fix_pt)(nc))
				b.back() -= m_all[nc.l+nr.l*lu] * (this->*fix_vl)(nc);

			if (b.back() != 0) adjust=true;
	}

	// insert adjustments to non-fixed values
	int nu(iu_nf.size()),*iu(iu_nf.data());
	dcomplex *d(b.data());
	if (adjust) VecSetValues(Bu,nu,iu,d,ADD_VALUES);
}

/**
 * function for calculating and inserting portions
 * of pressure source vector stemming from the zeroing out
 * of columns corresponding to fixed Dirichlet boundary conditions
 *
 * @param[in] eu velocity element in which to get matrix
 * @param[in] ep pressure element in which to get matrix
 * @param[in] da_pre pressure node DMDA
 * @param[out] Bp pressure source vector
 */
void brinkman::add_Diri_Bp(element &eu,element &ep,DM da_pre,Vec Bp) {

	// non-fixed mat and source values
	std::vector<dcomplex> b(lp);
	for(int i=0;i<lp;i++) b[i]=0;

	// get matrix values (row vel, col pre)
	dcomplex *mx(new dcomplex[2*lu*lp]),*my(mx+lu*lp);
	get_Auvp(eu,mx,my);

	// get global indices
	std::vector<int> gip(lp);
	get_p_global_indices(ep,da_pre,gip.data());

	// loop through pressure (row) nodes, vel (col) nodes, storing if fixed
	bool adj = false;
	for (node nr:ep) {
		for (node nc:eu) {
			if (fix_pt_ux(nc)) b[nr.l] -= mx[nr.l+nc.l*lp] * fix_vl_ux(nc);
			if (fix_pt_uy(nc)) b[nr.l] -= my[nr.l+nc.l*lp] * fix_vl_uy(nc);
		}

		if (b[nr.l] != 0) adj=true;
	}

	// cleanup
	delete[] mx;

	// adjust source values if we have any nonzeros
	int np(gip.size()),*ip(gip.data());
	dcomplex *d(b.data());
	if (adj) VecSetValues(Bp,np,ip,d,ADD_VALUES);
}

/**
 * Calculate and insert portion of velocity source vector
 * based on body force struct
 *
 * @param[in] eu velocity element
 * @param[in] da_vel velocity node DMDA
 * @param[out] Bu velocity source vector
 */
void brinkman::add_BU(element &eu,DM da_vel,Vec Bu) {

	// grab data
	std::vector<dcomplex> fx(lu),fy(lu);
	get_BU(eu,fx.data(),fy.data());

	// pull in global indices
	std::vector<int> iu(lu),iv(lu);
	get_u_global_indices(eu,da_vel,iu.data());
	get_v_global_indices(eu,da_vel,iv.data());

	// set values
	VecSetValues(Bu,lu,iu.data(),fx.data(),ADD_VALUES);
	VecSetValues(Bu,lu,iv.data(),fy.data(),ADD_VALUES);
}

/**
 * Function for inserting values into the velocity source vector 
 * corresponding to fixed Dirichlet boundary conditions (i.e. in
 * empty rows with 1s on the diagonal)
 *
 * @param[in] eu the velocity element in which to fix matrix
 * @param[in] da_vel velocity node DMDA
 * @param[out] Bu velocity source vector
 * @param[in] scale optional magnitude of diagonal matrix elements
 */
void brinkman::fix_BU(element &eu,DM da_vel,Vec Bu,const dcomplex scale) {

	// grab local vel inds
	std::vector<int> iu(lu),iv(lu);
	get_u_local_indices(eu,da_vel,iu.data());
	get_v_local_indices(eu,da_vel,iv.data());

	// loop through nodes
	for(node n: eu) {

		// x fixes
		if (fix_pt_ux(n)) {

			// grab value, scaled value, and index
			dcomplex v(fix_vl_ux(n)),sc_v(scale*v);
			int loc=iu[n.l];

			// fix mat diagonal, solution and source
			VecSetValuesLocal(Bu,1,&loc,&sc_v,INSERT_VALUES);
		}

		// y fixes
		if (fix_pt_uy(n)) {

			// grab value, scaled value, and index
			dcomplex v(fix_vl_uy(n)),sc_v(scale*v);
			int loc=iv[n.l];

			// fix mat diagonal, solution and source
			VecSetValuesLocal(Bu,1,&loc,&sc_v,INSERT_VALUES);
		}
	}
}

/**
 * Function for placing values on the diagonal of the velocity-velocity
 * submatrix to reflect a Dirichlet boundary condition
 *
 * @param[in] eu the velocity element in which to fix matrix
 * @param[in] da_vel velocity node DMDA
 * @param[out] Auu local velocity-velocity submatrix
 * @param[in] scale optional magnitude of diagonal matrix elements
 */
void brinkman::fix_AUU(element &eu,DM da_vel,Mat Auu,const dcomplex scale) {

	// grab local vel inds
	std::vector<int> iu(lu),iv(lu);
	get_u_local_indices(eu,da_vel,iu.data());
	get_v_local_indices(eu,da_vel,iv.data());

	// loop through nodes
	for(node n: eu) {
		if (fix_pt_ux(n)) {

			// grab value, scaled value, and index
			dcomplex v(fix_vl_ux(n)),sc_v(scale*v);
			int loc=iu[n.l];

			// fix mat diagonal, solution and source
			MatSetValuesLocal(Auu,1,&loc,1,&loc,&scale,INSERT_VALUES);
		}

		if (fix_pt_uy(n)) {

			// grab value, scaled value, and index
			dcomplex v(fix_vl_uy(n)),sc_v(scale*v);
			int loc=iv[n.l];

			// fix mat diagonal, solution and source
			MatSetValuesLocal(Auu,1,&loc,1,&loc,&scale,INSERT_VALUES);
		}
	}
}

/**
 * Helper function to identify the level of coarsening relative to
 * this brinkman given a set of velocity and pressure grid sizes
 *
 * @param[in] Mu velocity nodes in x-direction
 * @param[in] Nu velocity nodes in y-direction
 * @param[in] Mp pressure nodes in x-direction
 * @param[in] Np pressure nodes in y-direction
 * @param[out] cx level of coarsening in the x-direction
 * @param[out] cy level of coarsening in the y-direction
 */
void brinkman::coarsen_factors_vel(int Mu,int Nu,int &cx,int &cy) {

	// fine grid sizes
	int fmu=ump,fnu=unp;

	// coarse grid sizes
	int cmu=Mu, cnu=Nu;

	// adjust for periodicity
	if (!hp) {fmu--; cmu--;}
	if (!vp) {fnu--; cnu--;}

	// check for exact multiple
	if (fmu%cmu!=0) bail(-1,"x vel not an exact multiple for coarsening "
		"(%d -> %d), (%d -> %d) adjusted for periodicity\n",ump,Mu,fmu,cmu);
	if (fnu%cnu!=0) bail(-1,"y vel not an exact multiple for coarsening "
		"(%d -> %d), (%d -> %d) adjusted for periodicity\n",unp,Nu,fnu,cnu);

	// set factors
	cx = fmu/cmu; cy = fnu/cnu;
}

/**
 * Helper function to identify the level of coarsening relative to
 * this brinkman given a set of velocity and pressure grid sizes
 *
 * @param[in] Mu velocity nodes in x-direction
 * @param[in] Nu velocity nodes in y-direction
 * @param[in] Mp pressure nodes in x-direction
 * @param[in] Np pressure nodes in y-direction
 * @param[out] cx level of coarsening in the x-direction
 * @param[out] cy level of coarsening in the y-direction
 */
void brinkman::coarsen_factors_pre(int Mp,int Np,int &cx,int &cy) {

	// fine grid sizes
	int fmp=pmp,fnp=pnp;

	// coarse grid sizes
	int cmp=Mp, cnp=Np;

	// adjust for periodicity
	if (!hp) {fmp--; cmp--;}
	if (!vp) {fnp--; cnp--;}

	// check for exact multiple
	if (fmp%cmp!=0) bail(-1,"x pre not an exact multiple for coarsening "
		"(%d -> %d), (%d -> %d) adjusted for periodicity\n",pmp,Mp,fmp,cmp);
	if (fnp%cnp!=0) bail(-1,"y pre not an exact multiple for coarsening "
		"(%d -> %d), (%d -> %d) adjusted for periodicity\n",pnp,Np,fnp,cnp);

	// set factors
	cx = fmp/cmp; cy = fnp/cnp;
}

/**
 * Helper function to identify the level of coarsening relative to
 * this brinkman given a set of velocity and pressure grid sizes
 *
 * @param[in] Mu velocity nodes in x-direction
 * @param[in] Nu velocity nodes in y-direction
 * @param[in] Mp pressure nodes in x-direction
 * @param[in] Np pressure nodes in y-direction
 * @param[out] cx level of coarsening in the x-direction
 * @param[out] cy level of coarsening in the y-direction
 */
void brinkman::coarsen_factors(int Mu,int Nu,
		int Mp,int Np,int &cx,int &cy) {

	// fine grid sizes
	int fmu=ump,fnu=unp,fmp=pmp,fnp=pnp;

	// coarse grid sizes
	int cmu=Mu, cnu=Nu, cmp=Mp, cnp=Np;

	// adjust for periodicity
	if (!hp) {fmu--; fmp--; cmu--; cmp--;}
	if (!vp) {fnu--; fnp--; cnu--; cnp--;}

	// check for exact multiple
	if (fmu%cmu!=0) bail(-1,"x vel not an exact multiple for coarsening "
		"(%d -> %d), (%d -> %d) adjusted for periodicity\n",ump,Mu,fmu,cmu);
	if (fnu%cnu!=0) bail(-1,"y vel not an exact multiple for coarsening "
		"(%d -> %d), (%d -> %d) adjusted for periodicity\n",unp,Nu,fnu,cnu);
	if (fmp%cmp!=0) bail(-1,"x pre not an exact multiple for coarsening "
		"(%d -> %d), (%d -> %d) adjusted for periodicity\n",pmp,Mp,fmp,cmp);
	if (fnp%cnp!=0) bail(-1,"y pre not an exact multiple for coarsening "
		"(%d -> %d), (%d -> %d) adjusted for periodicity\n",pnp,Np,fnp,cnp);

	// set factors
	cx = fmu/cmu; cy = fnu/cnu;

	// check vel/pre have same factors
	if (fmp/cmp != cx) bail(-1,"x vel and pre factors inconsistent "
		"(%d vel, %d pre)\n",cx,fmp/cmp);
	if (fnp/cnp != cy) bail(-1,"y vel and pre factors inconsistent "
		"(%d vel, %d pre)\n",cy,fnp/cnp);
}

/**
 * Given the provided DM composite, allocates and returns a
 * brinkman corresponding to the coarsened grid
 *
 * @param[in] comp DMComposite of velocity and pressure grids
 * @return a brinkman for the coarsened grids
 */
brinkman *brinkman::coarsen(DM comp) {

	// grid sizes
	int cx,cy;
	DM vel,pre;
	DMDALocalInfo vi,pi;
	switch (dm_type(comp)) {

		case dmgr::sys:

			// grab individual vel and pres grids
			DMCompositeGetEntries(comp,&vel,&pre);

			// grab infos
			DMDAGetLocalInfo(vel,&vi);
			DMDAGetLocalInfo(pre,&pi);

			coarsen_factors(vi.mx,vi.my,pi.mx,pi.my,cx,cy);
			break;

		case dmgr::vel:

			vel=comp;
			DMDAGetLocalInfo(vel,&vi);
			coarsen_factors_vel(vi.mx,vi.my,cx,cy);
			break;

		case dmgr::pre:

			pre=comp;
			DMDAGetLocalInfo(pre,&pi);
			coarsen_factors_pre(pi.mx,pi.my,cx,cy);
			break;

		default:

			bail(0,"bad DM for coarsening function\n");
	}


	// call virtual func based on grid sizes
	return coarsen(cx,cy);
}

///////////////// brinkman info (driven_cyl) ////////////////////

/**
 * Allocates a new info class coarsenedd with the provided factors
 * in each direction
 *
 * @param[in] cx
 * @param[in] cy
 */
brinkman *driven_cyl::coarsen(int cx,int cy) {

	// allocate new info, don't delete source func
	driven_cyl *par = new driven_cyl(po,m/cx,n/cy,hp,vp,h,L,f,mu,alpha,V_inf);
	par->del_src=false;

	return par;
}

///////////////// brinkman info (one_sphere) ////////////////////

/**
 * base constructor
 *
 * @param[in] p poly order
 * @param[in] m number of horizontal elements
 * @param[in] n number of vertical elements
 * @param[in] a inner sphere radius
 * @param[in] R outer sphere radius
 * @param[in] f body-force source function
 * @param[in] mu viscosity
 * @param[in] alph inverse drag length
 * @param[in] V_i inner sphere vel
 * @param[in] V_o outer sphere vel
 */
one_sphere::one_sphere(int p,int m,int n,double a_,double R_,
	tsrc_base<vcomplex> *f,double mu,double alph,vcomplex Vi,vcomplex Vo):
	brinkman(ptype::one_sphere,
		p,m,n,false,true,new spherical(m,n,a_,R_),f,mu,alph,Vo),del_src(true),
	a(a_),R(R_),V_i(Vi),dV_i(new vcomplex[unp]),dV_o(new vcomplex[unp]) {

	// zero all variable density
	for(int i=0;i<unp;i++) dV_i[i]=dV_o[i]=vcomplex::zero;
}

/**
 * steady solution induced by leading-order brinkman solution
 *
 * @param[in] par one_sphere problem parameters
 * @param[in] br oscillatory Brinkman solution
 * @param[in] type which parts of steady solution to take
 *   (0: whole thing, 1: only boundary vel, 2: only Reynolds stress)
 */
one_sphere::one_sphere(one_sphere *par,int type):
	one_sphere(par->po,par->m,par->n,par->a,par->R,
		type==1 ?
		static_cast<tsrc_base<vcomplex>*>(new homog()) :
		static_cast<tsrc_base<vcomplex>*>(new reynolds(par)),
	par->mu,0,vcomplex::zero,vcomplex::zero) {

	if (type==2) return;

	// loop through boundary nodes
//	vcomplex v(vcomplex::zero);
	for(int gj=0,gi=0;gj<unp;gj++) {

		// grab node, assign boundary vel
		node n(0,0,0,gi+gj*ump,gi,gj);
		dV_i[n.gj] = par->steady_bvel(n);
	}

}

/**
 * class destructor
 */
one_sphere::~one_sphere() {
	if (del_src) delete f;
	delete dm;
	delete[] dV_i;
	delete[] dV_o;
}

/**
 * save the one_sphere struct as a binary
 *
 * @param[in] fn filename
 * @return number of bytes written
 */
void one_sphere::save_binary(const char *fn) {

	// open writefile
	FILE *fb = fopen(fn,"wb");

	// write sim parameters
	fwrite(&type,    sizeof(ptype),  1,fb);	
	fwrite(&po,      sizeof(int),    1,fb);
	fwrite(&m,       sizeof(int),    1,fb);
	fwrite(&n,       sizeof(int),    1,fb);
	fwrite(&a,       sizeof(double), 1,fb);
	fwrite(&R,       sizeof(double), 1,fb);
	fwrite(&mu,      sizeof(double), 1,fb);
	fwrite(&alpha,   sizeof(double), 1,fb);

	// write boundary conditions
	fwrite(&V_i,   sizeof(vcomplex),1,fb);
	fwrite(&V_inf, sizeof(vcomplex),1,fb);
	fwrite(dV_i,   sizeof(vcomplex),unp,fb);
	fwrite(dV_o,   sizeof(vcomplex),unp,fb);

	fwrite(&data,  sizeof(bool),1,fb);
	if (data) write_data(fb);

	fclose(fb);
}

/**
 * load a one_sphere struct from a binary
 *
 * @param[in] fn filename
 * @param[out] number of bytes in one_sphere binary
 * @return newly allocated one_sphere struct
 */
one_sphere* one_sphere::load_binary(const char *fn) {

	// open up a binary file
	FILE *fb = fopen(fn,"rb");

	// set up variables to read in
	ptype pt;
	// bc_type bc;
	int p,m,n;
	double a,R,mu,alpha;

	// check if this is the right type
	fread(&pt,    sizeof(ptype), 1,fb);
	assert(pt==ptype::one_sphere);

	// read in problem information
	fread(&p,     sizeof(int),    1,fb);
	fread(&m,     sizeof(int),    1,fb);
	fread(&n,     sizeof(int),    1,fb);
	fread(&a,     sizeof(double), 1,fb);
	fread(&R,     sizeof(double), 1,fb);
	fread(&mu,    sizeof(double), 1,fb);
	fread(&alpha, sizeof(double), 1,fb);

	// allocate a new object
	// use const vel as default, then read in boundary vals
	one_sphere *par = new one_sphere(p,m,n,a,R,mu,alpha,1.);

	// read in boundary values
	fread(&par->V_i,   sizeof(vcomplex),1,fb);
	fread(&par->V_inf, sizeof(vcomplex),1,fb);
	fread(par->dV_i,   sizeof(vcomplex),par->unp,fb);
	fread(par->dV_o,   sizeof(vcomplex),par->unp,fb);

	bool has_dat;
	fread(&has_dat,  sizeof(bool),1,fb);
	if (has_dat) par->read_data(fb);

	// store the spot in the file, close
	fclose(fb);

	return par;
}

/**
 * returns the boundary velocity induced by steady streaming
 *
 * @param[in] nb boundary node at which to calculate boundary vel
 * @param[in] br brinkman solution with oscillatory solution
 * @return the steady velocity at the node from the brinkman solution
 */
vcomplex one_sphere::steady_bvel(node &nb) {

	// only do this on the inner sphere
	assert(nb.gi==0);

	// element index depends on position in global grid
	int ei = 0;
	int ej = nb.gj / po;

	// we start by assume we're on low end of an element

	// node x-index
	int ni = 0;
	int nj = nb.gj % po;
	int nn = ni + nj*(po+1);

	double Ux,Uy;
	vcomplex retval(0);

	// grab referent to element / node
	element e = gr_u(ei,ej);
	node nd(nn,ni,nj);

	// get velocity there
	vcomplex uB = (0.5*I*u_grad_u(e,nd)).real();
	retval = uB;

	// we're done if this is not a boundary
	if (nj!=0) return retval;

	// otherwise, average with next one
	retval = 0.5 * retval;

	// grab next element
	ej--; if (ej<0) ej=n-1;
	nj = po;
	nn = ni + nj*(po+1);

	// grab new element and node
	e = gr_u(ei,ej);
	nd = node(nn,ni,nj);

	// grab new val, average it in
	uB = (0.5*I*u_grad_u(e,nd)).real();
	retval += 0.5*uB;

	return retval;
}

/**
 * Return a new brinkman coarsened in each direction
 * according to the provided factor in each direction
 *
 * @param[in] cx factor to coarsen in x-direction
 * @param[in] cy factor to coarsen in y-direction
 */
brinkman* one_sphere::coarsen(int cx,int cy) {

	// allocate new one_sphere, but don't delete source func
	one_sphere *par = new one_sphere(po,m/cx,n/cy,a,R,f,mu,alpha,V_i,V_inf);
	par->del_src=false;

	// coarsen the boundary vels
	for(int i=0; i<par->unp; i++) {
		par->dV_i[i] = dV_i[i*cy];
		par->dV_o[i] = dV_o[i*cy];
	}

	return par;
}

///////////////// brinkman info (two_spheres) ////////////////////

/**
 * base constructor
 *
 * @param[in] p poly order
 * @param[in] m number of horizontal elements
 * @param[in] n number of vertical elements
 * @param[in] hp horizontal periodicity
 * @param[in] vp vertical periodicity
 * @param[in] D sphere separation
 * @param[in] r1 left sphere radius
 * @param[in] r2 right sphere radius
 * @param[in] f body-force source function
 * @param[in] mu viscosity
 * @param[in] alph inverse drag length
 * @param[in] Vx_l left sphere x vel
 * @param[in] Vy_l left sphere y vel
 * @param[in] Vx_r right sphere x vel
 * @param[in] Vy_r right sphere y vel
 */
two_spheres::two_spheres(int p,int m,int n,
	double D,double r1,double r2,tsrc_base<vcomplex> *f,double mu,double alph,
	vcomplex V_l_,vcomplex V_r_,vcomplex V_inf):
	brinkman(ptype::two_spheres,p,m,n,false,true,

		new bispherical(

		
			m/4,(m/2)-(m/4), (m/2)-(m/4), m/4,
			n,D,r1,r2),//,alph),

		f,mu,alph,V_inf),
	del_src(true),rl(r1),rr(r2),d(D),
	xl(bipolar_half_plane::drr2a(D,r1,r2)
			/tanh(bipolar_half_plane::drr2xi1(D,r1,r2))),
	xr(bipolar_half_plane::drr2a(D,r1,r2)
			/tanh(bipolar_half_plane::drr2xi2(D,r1,r2))),
	xil(bipolar_half_plane::drr2xi1(D,r1,r2)),
	xir(bipolar_half_plane::drr2xi2(D,r1,r2)),
	V_l(V_l_),V_r(V_r_),
	dV_l(new vcomplex[unp]),dV_r(new vcomplex[unp]),F_a(vcomplex::zero) {

	// initialize all boundary vels to zero
	for(int i=0;i<unp;i++) dV_l[i]=dV_r[i]=vcomplex::zero;
}

/**
 * steady solution induced by leading-order brinkman solution
 *
 * @param[in] par two_sphere problem parameters
 * @param[in] br oscillatory Brinkman solution
 * @param[in] V swimming velocity
 * @param[in] type which parts of steady solution to take
 *   (0: whole thing, 1: only boundary vel, 2: only Reynolds stress, 
 *    3: only net osc force)
 */
two_spheres::two_spheres(two_spheres* par,dcomplex V,int type):
	two_spheres(par->po,par->m,par->n,par->d,par->rl,par->rr,
		(type==1 || type==3) ?
		static_cast<tsrc_base<vcomplex>*>(new homog()) :
		static_cast<tsrc_base<vcomplex>*>(new reynolds(par)),
	par->mu,0,vcomplex::zero,vcomplex::zero,vcomplex(-V,0)) {

	// don't generate boundary velocities or net stress
	// if looking at Reynolds stress-only
	if (type==0 || type==1){ 

		// loop through boundary nodes
		for(int gj=0,gi_l=0,gi_r=(ump-1);gj<unp;gj++) {

			// grab left and right-sphere boundary nodes
			node nl(0,0,0,gi_l+gj*ump,gi_l,gj);
			node nr(0,0,0,gi_r+gj*ump,gi_r,gj);

			// assign boundary velocities
			dV_l[nl.gj] = par->steady_bvel(nl);
			dV_r[nr.gj] = par->steady_bvel(nr);
		}
	}
}

/**
 * returns the boundary velocity induced by steady streaming
 *
 * @param[in] nb boundary node at which to calculate boundary vel
 * @param[in] br brinkman solution with oscillatory solution
 * @return the steady velocity at the node from the brinkman solution
 */
vcomplex two_spheres::steady_bvel(node &nb) {

	// only do this on the two spheres
	assert(nb.gi==0 || nb.gi==(ump-1));

	// element index depends on position in global grid
	int ei = nb.gi==0 ? 0 : (m-1);
	int ej = nb.gj / po;

	// we start by assume we're on low end of an element

	// node x-,y-,local index
	int ni = nb.gi==0 ? 0 : po;
	int nj = nb.gj % po;
	int nn = ni + nj*(po+1);

	double Ux,Uy;
	vcomplex retval(0);

	// grab referent to element / node
	element e = gr_u(ei,ej);
	node nd(nn,ni,nj);

	// get velocity there
	vcomplex uB = (0.5*I*u_grad_u(e,nd)).real();
	retval = uB;//vcomplex(Ux,Uy);

	// we're done if this is not an element boundary
	if (nj!=0) return retval;

	// otherwise, average with next one
	retval = 0.5 * retval;

	// grab next element
	ej--; if (ej<0) ej=n-1;
	nj = po;
	nn = ni + nj*(po+1);

	// grab new element and node
	e = gr_u(ei,ej);
	nd = node(nn,ni,nj);

	// grab new val, average it in
	uB = (0.5*I*u_grad_u(e,nd)).real();
	retval += 0.5*uB;

	return retval;
}

/*
 * class destructor
 */
two_spheres::~two_spheres() {
	if (del_src) delete f;
	delete dm;
	delete[] dV_l;
	delete[] dV_r;
}

/**
 * save the two_spheres struct as a binary
 *
 * @param[in] fn filename
 * @return number of bytes written
 */
void two_spheres::save_binary(const char *fn) {

	// open writefile
	FILE *fb = fopen(fn,"wb");

	// write simulation parameters
	fwrite(&type,    sizeof(ptype),  1,fb);	
	fwrite(&po,      sizeof(int),    1,fb);
	fwrite(&m,       sizeof(int),    1,fb);
	fwrite(&n,       sizeof(int),    1,fb);
	fwrite(&d,       sizeof(double), 1,fb);
	fwrite(&rl,      sizeof(double), 1,fb);
	fwrite(&rr,      sizeof(double), 1,fb);
	fwrite(&mu,      sizeof(double), 1,fb);
	fwrite(&alpha,   sizeof(double), 1,fb);

	// write boundary vel data
	fwrite(&V_l,    sizeof(vcomplex),1,fb);
	fwrite(&V_r,    sizeof(vcomplex),1,fb);
	fwrite(&V_inf,  sizeof(vcomplex),1,fb);
	fwrite(dV_l,    sizeof(vcomplex),unp,fb);
	fwrite(dV_r,    sizeof(vcomplex),unp,fb);

	fwrite(&data,   sizeof(bool),1,fb);
	if (data) write_data(fb);

	fclose(fb);
}

/**
 * load a two_spheres struct from a binary
 *
 * @param[in] fn filename
 * @param[out] number of bytes in one_sphere binary
 * @return newly allocated two_spheres struct
 */
two_spheres* two_spheres::load_binary(const char *fn) {

	// open up a binary file
	FILE *fb = fopen(fn,"rb");

	// set up variables to read in
	ptype pt;
	// bc_type bc;
	int p,m,n;
	double d,rl,rr,mu,alpha;

	// check if this is the right type
	fread(&pt,    sizeof(ptype), 1,fb);
	assert(pt==ptype::two_spheres);

	// read system parameters
	fread(&p,     sizeof(int),    1,fb);
	fread(&m,     sizeof(int),    1,fb);
	fread(&n,     sizeof(int),    1,fb);
	fread(&d,     sizeof(double), 1,fb);
	fread(&rl,    sizeof(double), 1,fb);
	fread(&rr,    sizeof(double), 1,fb);
	fread(&mu,    sizeof(double), 1,fb);
	fread(&alpha, sizeof(double), 1,fb);

	// allocate
	two_spheres *par = new two_spheres(p,m,n,d,rl,rr,mu,alpha,0,0);

	// read boudnary data
	fread(&par->V_l,   sizeof(vcomplex),1,fb);
	fread(&par->V_r,   sizeof(vcomplex),1,fb);
	fread(&par->V_inf, sizeof(vcomplex),1,fb);
	fread(par->dV_l,   sizeof(vcomplex),par->unp,fb);
	fread(par->dV_r,   sizeof(vcomplex),par->unp,fb);

	bool has_dat;
	fread(&has_dat,    sizeof(bool),1,fb);
	if (has_dat) par->read_data(fb);

	// store the spot in the file, close
	fclose(fb);

	// allocate a new object
	return par;
}

/**
 * given a brinkman solution, calculate the viscous drag
 * force applied to one of the spheres by the fluid
 *
 * @param[in] br brinkman solution
 * @param[in] left whether the left sphere (right if false)
 */
dcomplex two_spheres::drag_force(bool left,bool approx) {

	if (approx) {
		dcomplex V = left?V_l.x:V_r.x;
		dcomplex v = left?V_r.x:V_l.x;

		double R = left?xl:xr;
		double r = left?xr:xl;

		V -= V_inf.x;
		v -= V_inf.x;
		v *= (1.5 * (r/d));

		dcomplex F =-6 * M_PI * R * (V-v);
		
		return F;
	}

	// total up infinitesimal forces
	dcomplex F(0);

	// 1D gauss rule
	gauss_quad_1d q((gr_u.p+1)/2 + 3);

	// coord trafos
	trafo_base *ct_p,*ct_u;

	// u,p values and derivatives
	dcomplex uu,uu_x,uu_y,vv,vv_x,vv_y,pp,pp_x,pp_y;

	// solution datapoints
	dcomplex *usol = new dcomplex[lu],
		*vsol = new dcomplex[lu],*psol = new dcomplex[lp];

	// Jacobean info / normal vector
	double det;
	double x,y;
	double xi_x,et_x,xi_y,et_y,x_et,y_et;
	double nx,ny,n2;

	// loop over xi_1 and xi_2, i.e. left/right sides of cartesian grid
	for (int j=0;j<gr_u.n;j++) {

		// element index
		int i = left ? 0 : (gr_u.m-1);

		// left/right elements
		element eu = gr_u(i,j);
		element ep = gr_p(i,j);

		// trafos
		ct_u = dm->alloc_element_geometry(eu);
		ct_p = dm->alloc_element_geometry(ep);

		// grab vel solution
		for (node nu:eu) {
			usol[nu.l] = ux[nu.g];
			vsol[nu.l] = uy[nu.g];
		}

		// ...and p solution
		for (node np:ep) psol[np.l] = p[np.g];

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

			// construct line segment length
			x_et =-xi_y/tmp_det;
			y_et = xi_x/tmp_det; 

			// relevant determinant is 1D in xi-direction
			det = sqrt(x_et*x_et + y_et*y_et) * fabs(y) * M_PI;
			
			// construct normal unit vector
			// from xi-gradient
			nx = xi_x; ny = xi_y;
			if (!left) {nx *= -1; ny *= -1;}
			
			// normalize
			n2 = nx*nx + ny*ny;
			nx /= sqrt(n2);
			ny /= sqrt(n2);

			// pull out vel and pressure sol'ns
			mast_u.eval<dcomplex>(*ct_u,xi,et,usol,uu,uu_x,uu_y);
			mast_u.eval<dcomplex>(*ct_u,xi,et,vsol,vv,vv_x,vv_y);
			mast_p.eval<dcomplex>(*ct_p,xi,et,psol,pp,pp_x,pp_y);

			// construct stress tensor
			dcomplex s_xx = 2*mu* uu_x - pp;
			dcomplex s_xy =   mu*(uu_y + vv_x);

			// force is x-component of stress
			F += (s_xx * nx + s_xy * ny)*q.w[g]*det;
		}

		dm->free_element_geometry(ct_u);
		dm->free_element_geometry(ct_p);
	}

	delete[] psol;
	delete[] vsol;
	delete[] usol;
	
	return F;

}

/**
 * given a brinkman solution, calculate the viscous drag
 * force applied to one of the spheres by the fluid
 *
 * @param[in] br brinkman solution
 * @param[in] left whether the left sphere (right if false)
 *
 */
dcomplex two_spheres::avg_vel(bool left) {


	// total up x-velocity
	dcomplex V(0);
	double area(0);

	// 1D gauss rule
	gauss_quad_1d q((gr_u.p+1)/2 + 3);

	// coord trafos
	trafo_base *ct_u;

	// solution datapoints
	dcomplex *usol = new dcomplex[lu];
	dcomplex *vsol = new dcomplex[lu];

	dcomplex uu,vv;

	// Jacobean info / normal vector
	double det;
	double x,y;
	double xi_x,et_x,xi_y,et_y,x_et,y_et;

	// radius
	double R = left?rl:rr;

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

			if ( ( left && (nu.li != 0) ) || ( !left && (nu.li != po)) ) {
				usol[nu.l] = 0; 
				vsol[nu.l] = 0;
			} else  {
				usol[nu.l] = left ? dV_l[nu.gj].x : dV_r[nu.gj].x;
				vsol[nu.l] = left ? dV_l[nu.gj].y : dV_r[nu.gj].y;
			}
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

			// normal vector
			double nx = x-(left?xl:xr);
			double ny = y;

			double th = atan2(ny,nx);

			double n_mag = sqrt(nx*nx + ny*ny);
			nx /= n_mag; ny /= n_mag;

			// construct line segment length
			x_et =-xi_y/tmp_det;
			y_et = xi_x/tmp_det; 

			// relevant determinant is 1D in xi-direction
			det = sqrt(x_et*x_et + y_et*y_et) * fabs(y) * M_PI;
			
			// pull out vel and pressure sol'ns
			dcomplex _1,_2;
			mast_u.eval<dcomplex>(*ct_u,xi,et,usol,uu,_1,_2);
			mast_u.eval<dcomplex>(*ct_u,xi,et,vsol,vv,_1,_2);

			dcomplex u_r = uu*nx + vv*ny;

			V += uu * q.w[g] * det;
			area += q.w[g] * det;
		}

		dm->free_element_geometry(ct_u);
	}

	delete[] usol;


	V /= area;
	
	return V;
}

/**
 * Perform a series of solves on the flow field such that the
 * applied forces on each sphere are equal and opposite using
 * Broyden's method
 *
 * @param[in] br the brinkman solution subject to the solve
 * @param[in] Fmag the magnitude of the equal/opposite forces
 */
void two_spheres::solve_force_balance(dcomplex Fmag) {

	petsc_solve petsc(this,sol_type);

	// set max iterations, solution tolerance
	static const int max_its = 100;
	static const double tol = std::numeric_limits<double>::epsilon()*1e4;
	
	// previous/new velocity, previous force values, differences
	dcomplex Vl_old,Vr_old,Fl_old,Fr_old,Vl_new,Vr_new,dVl,dVr,dFl,dFr;
	
	// grab velocity references, force vals
	dcomplex &Vl=V_l.x,&Vr=V_r.x,Fl,Fr;

	// do 2 solves to estimate Jacobean


	// try changing left velocity, get Jacobean from resulting
	// force values (V=0 => F=0 so get third point for free,
	// use magnitude 1 so don't need to divide by Vl)
	// derivs: col index so dF = J.dV
	Vl = 1; Vr = 0;
	petsc.solve();
	dcomplex Jll(appl_force(true)),Jrl(appl_force(false));
	

	// try changing right velocity
	Vl = 0; Vr = 1;
	petsc.solve();
	dcomplex Jlr(appl_force(true)),Jrr(appl_force(false));

	// invert Jacobean
	dcomplex Jd(Jll*Jrr-Jlr*Jrl);
	dcomplex Jill(Jrr/Jd),Jilr(-Jlr/Jd),Jirl(-Jrl/Jd),Jirr(Jll/Jd);
	
	// use inverse Jacobean to start serach at the
	// vel corresponding to Fl = -F_mag, Fr = F_mag
	Vl_new = -Jill*Fmag+Jilr*Fmag;
	Vr_new = -Jirl*Fmag+Jirr*Fmag;
	Vl=Fl=Vr=Fr=0;

	int its=0; //count iterations
	do {

		// store old values
		Vl_old=Vl; Vr_old=Vr; Fl_old=Fl; Fr_old=Fr;

		// full update
		Vl=Vl_new; Vr=Vr_new;

		petsc.solve();
		Fl = appl_force(true)+Fmag; Fr = appl_force(false)-Fmag;

		// record diffs
		dVl = Vl-Vl_old; dVr = Vr-Vr_old;
		dFl = Fl-Fl_old; dFr = Fr-Fr_old;

		// update to Jacobean: "good" Broyden rule

		// vector: dD = (dX - J^-1.df)/(dX.J^-1.df)
		dcomplex dUl(Jill*dFl+Jilr*dFr),dUr(Jirl*dFl+Jirr*dFr);
		dcomplex dot(dVl*dUl+dVr*dUr);
		dcomplex dDl((dVl-dUl)/dot),dDr((dVr-dUr)/dot);

		// update Jacobean
		Jill += dDl*dVl*Jill + dDl*dVr*Jirl;
		Jilr += dDl*dVl*Jilr + dDl*dVr*Jirr;
		Jirl += dDr*dVl*Jill + dDr*dVr*Jirl;
		Jirr += dDr*dVl*Jilr + dDr*dVr*Jirr;

		// new update
		Vl_new = Vl - Jill*Fl - Jilr*Fr;
		Vr_new = Vr - Jirl*Fl - Jirr*Fr;


	} while (std::norm(Fl)+std::norm(Fr) > 4*tol*tol && ++its<max_its);
}

/**
 * Perform a series of solves on the flow field such that there
 * is no net drag force on the two-sphere system 
 *
 * @param[in] br the brinkman solution subject to the solve
 */
void two_spheres::solve_force_free() {

	petsc_solve petsc(this,sol_type);

	// set max iterations
	static const int max_its = 100;
	static const double tol = std::numeric_limits<double>::epsilon()*1e4;

	// new and old, diffs, differences
	dcomplex V_new,V_old,F_old,dV,dF;

	// grab reference to velocity, do initial solve w/0
	dcomplex &V = V_inf.x; V = 0;
	petsc.solve();
	dcomplex F = drag_force(true)+drag_force(false) - F_a.x;


	// do another w/1
	V = 1;
	petsc.solve();


	dF = drag_force(true)+drag_force(false) - F_a.x - F;

	// set up first guess where we expect F=0
	V_new = -F/dF;

	int its=0; // count iterations
	do {

		// store old, do solve with new value
		V_old=V; F_old=F; V=V_new;
		petsc.solve();
		F = drag_force(true)+drag_force(false) - F_a.x;

		// store diffs
		dV=V-V_old; dF=F-F_old;

		// secant update
		V_new = V-F*dV/dF;

	} while (std::norm(F) > 2*tol*tol && ++its<max_its);
}

/**
 * Return a new brinkman coarsened in each direction
 * according to the provided factor in each direction
 *
 * @param[in] cx factor to coarsen in x-direction
 * @param[in] cy factor to coarsen in y-direction
 */
brinkman* two_spheres::coarsen(int cx,int cy) {

	// allocate new one_sphere
	two_spheres *par = new two_spheres(po,m/cx,n/cy,
			d,rl,rr,f,mu,alpha,V_l,V_r,V_inf);
	par->del_src = false;

	// coarsen the boundary vels
	for(int i=0; i<par->unp; i++) {
		par->dV_l[i] = dV_l[i*cy];
		par->dV_r[i] = dV_r[i*cy];
	}

	return par;
}

//////////////////////// petsc_solve //////////////////////////

/**
 * Create and allocate a PETSc solver
 */
petsc_solve::petsc_solve(brinkman *prob,method type) {

	assign(prob);

	// allocate the grid
	DM sys;
	alloc_sys_dmda(prob,sys);

	// create the solver
	KSPCreate(PETSC_COMM_WORLD,&ksp);
	KSPSetDM(ksp,sys);

	// note that only LU works
	switch (type) {
		case method::lu:    set_lu();    break;
	}

	// give functions to compute operators
	KSPSetComputeOperators(ksp,&calculate_A,(void*)&ctx);

	KSPSetFromOptions(ksp);
	KSPSetUp(ksp);
}

/**
 * Solve problem
 */
void petsc_solve::solve() {

	KSPSetComputeRHS(ksp,&calculate_B,(void*)&ctx);

	KSPSolve(ksp,NULL,NULL);

	check_sol();

	ctx.prob->alloc_data();
	ctx.prob->read(*this);
}

/**
 * Helper function to count nodes on processors
 */
void petsc_solve::proc_counts(brinkman *prob,
	std::vector<int> &lux,std::vector<int> &luy,
	std::vector<int> &lpx,std::vector<int> &lpy) {

	// processor decomp
	int mp,np;
	proc_decomp(prob,mp,np);

	// reserve spots in vectors
	lux.resize(mp); lpx.resize(mp);
	luy.resize(np); lpy.resize(np);

	// check same lengths
	if (mp != lpx.size() || np != lpy.size())
		bail(-1,"mismatched node count vectors\n");

	// element numbers, poly order, periodicity
	int me(prob->m),ne(prob->n),u_po(prob->po),p_po(u_po-1);
	bool xprd(prob->hp),yprd(prob->vp);

	// count x procs
	for(int ip=0;ip<mp;ip++) {

		// how many elements on the proc?
		int n_els((ip+1)*me/mp - ip*me/mp);

		// convert with polynomial order
		lux[ip] = n_els*u_po;
		lpx[ip] = n_els*p_po;

		// add node to the end if non-periodic
		if (!xprd && ip==mp-1) {lux[ip]++; lpx[ip]++;}
	}

	// count y procs
	for(int ip=0;ip<np;ip++) {

		// how many elements on the proc?
		int n_els((ip+1)*ne/np - ip*ne/np);

		// convert with polynomial order
		luy[ip] = n_els*u_po;
		lpy[ip] = n_els*p_po;

		// add node to the end if non-periodic
		if (!yprd && ip==np-1) {luy[ip]++; lpy[ip]++;}
	}
}

/**
 * Helper function to find the best processor decomposition
 */
void petsc_solve::proc_decomp(brinkman *prob,int &mp,int &np) {

	int size;
	MPI_Comm_size(PETSC_COMM_WORLD,&size);

	// check each mp/np combo such that mp*np=size, pick the combo
	// with the smallest amount of communication
	int mp_opt,np_opt,comm_min(prob->m*prob->n);
	mp = size+1;
	while (--mp >= 1) if (size%mp==0) {

		// grab y procs
		np = size/mp;

		// measure comm area
		int xcomms=mp; if(!prob->hp) xcomms--;
		int ycomms=np; if(!prob->vp) ycomms--;

		// total communication
		int comm = (xcomms * prob->n) + (ycomms * prob->m);

		// store if good if (comm < comm_min) {
		comm_min=comm;
		mp_opt=mp; np_opt=np;
	}

	// use best
	mp=mp_opt; np=np_opt;
}
	
/**
 * Allocate the grid for the system
 */
void petsc_solve::alloc_sys_dmda(brinkman *prob,DM &sys) {

	// u vs p polynomial order
	int u_po(prob->po),p_po(u_po-1);

	// number of elements, dof, x-/y-periodicity
	int me(prob->m),ne(prob->n),u_dof(2),p_dof(1);
	bool xprd(prob->hp),yprd(prob->vp);

	// helpers, DMDA boundary types, MPI_Comms
	DMBoundaryType none(DM_BOUNDARY_NONE),prd(DM_BOUNDARY_PERIODIC);
	DMDAStencilType box(DMDA_STENCIL_BOX);
	MPI_Comm world(PETSC_COMM_WORLD);

	// generate node counts for each processor
	std::vector<int> lux,luy,lpx,lpy;
	proc_counts(prob,lux,luy,lpx,lpy);
	int Mu(sum(lux)),Nu(sum(luy)),Mp(sum(lpx)),Np(sum(lpy));
	int mp(lux.size()),np(luy.size());

	// DMs
	DM vel,pre;

	// velocity DMDA
	DMDACreate2d(world,xprd?prd:none,yprd?prd:none,box,
		Mu,Nu,mp,np,u_dof,u_po,lux.data(),luy.data(),&vel);
	DMSetOptionsPrefix(vel,"vel_");
	DMSetUp(vel);

	// pressure DMDA
	DMDACreate2d(world,xprd?prd:none,yprd?prd:none,box,
		Mp,Np,mp,np,p_dof,p_po,lpx.data(),lpy.data(),&pre);
	DMSetOptionsPrefix(pre,"pre_");
	DMSetUp(pre);

	// composite
	DMCompositeCreate(world,&sys);
	DMSetOptionsPrefix(sys,"sys_");
	DMCompositeAddDM(sys,vel);
	DMCompositeAddDM(sys,pre);
	DMDASetFieldName(vel,0,"u");
	DMDASetFieldName(vel,1,"v");
	DMDASetFieldName(pre,0,"p");
	DMSetUp(sys);
}

/**
 * provides command-line args to use LU factorization
 */
void petsc_solve::set_lu(const char* lu_type) {
	add_opt("-ksp_type","preonly");
	add_opt("-pc_type","lu");
	add_opt("-pc_factor_mat_solver_type",lu_type);
}

/**
 * Calculate the RHS of the system
 */
PetscErrorCode petsc_solve::calculate_B(KSP ksp,Vec B,void *ctx) {

	VecZeroEntries(B);


	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	// grab brinkman
	brinkman *par = static_cast<user_context*>(ctx)->prob;

	// get DM and system size
	DM sys,vel,pre;
	KSPGetDM(ksp,&sys);

	// get coarser brinkman
	brinkman *coarse = par->coarsen(sys);

	// pull out other DMDAs
	DMCompositeGetEntries(sys,&vel,&pre);

	ISLocalToGlobalMapping map_u,map_p;
	DMGetLocalToGlobalMapping(vel,&map_u);
	DMGetLocalToGlobalMapping(pre,&map_p);


	IS *iss;
	DMCompositeGetGlobalISs(sys,&iss);

	// split vec into vel and pressure
	Vec Bu,Bp;
	VecGetSubVector(B,iss[0],&Bu);
	VecGetSubVector(B,iss[1],&Bp);

	VecSetLocalToGlobalMapping(Bu,map_u);
	VecSetLocalToGlobalMapping(Bp,map_p);

	DMDALocalInfo inf;
	DMDAGetLocalInfo(pre,&inf);

	int Mu,Nu,Mp,Np;
	DMDAGetInfo(vel,NULL,
		&Mu,&Nu,NULL,
		NULL,NULL,NULL,
		NULL,NULL,NULL,
		NULL,NULL,NULL);
	DMDAGetInfo(pre,NULL,
		&Mp,&Np,NULL,
		NULL,NULL,NULL,
		NULL,NULL,NULL,
		NULL,NULL,NULL);
	

	// get element range
	int euai,eubi,euaj,eubj;
	int epai,epbi,epaj,epbj;
	brinkman::dm_fe_range(vel,par->po,  euai,eubi,euaj,eubj);
	brinkman::dm_fe_range(pre,par->po-1,epai,epbi,epaj,epbj);

	if (euai!=epai || eubi!=epbi || euaj!=epaj || eubj!=epbj) {
		bail(-1,"element range mismatch:\n"
				"P: [%d,%d) x [%d,%d)\n"
				"U: [%d,%d) x [%d,%d)\n",
				epai,epbi,epaj,epbj,euai,eubi,euaj,eubj);
	}


	for (int ej=euaj; ej<eubj; ej++) for (int ei=euai; ei<eubi; ei++) {

		element eu(coarse->gr_u(ei,ej)), ep(coarse->gr_p(ei,ej));

		coarse->add_BU(eu,vel,Bu);
		coarse->add_Diri_Bu(eu,vel,Bu);
		coarse->add_Diri_Bv(eu,vel,Bu);
		coarse->add_Diri_Bp(eu,ep,pre,Bp);
	}

	VecAssemblyBegin(Bu);
	VecAssemblyBegin(Bp);
	VecAssemblyEnd(Bp);
	VecAssemblyEnd(Bu);

	VecRestoreSubVector(B,iss[0],&Bu);
	VecRestoreSubVector(B,iss[1],&Bp);

	VecAssemblyBegin(B);
	VecAssemblyEnd(B);


	VecGetSubVector(B,iss[0],&Bu);
	VecSetLocalToGlobalMapping(Bu,map_u);

	for (int ej=euaj; ej<eubj; ej++) for (int ei=euai; ei<eubi; ei++) {

		element eu(coarse->gr_u(ei,ej)), ep(coarse->gr_p(ei,ej));
		coarse->fix_BU(eu,vel,Bu);
	}

	VecAssemblyBegin(Bu);
	VecAssemblyEnd(Bu);

	VecRestoreSubVector(B,iss[0],&Bu);


	VecAssemblyBegin(B);
	VecAssemblyEnd(B);



	return 0;
}

/**
 * Calculate the mass matrix of the pressure system
 */
PetscErrorCode petsc_solve::calculate_p_mass(KSP ksp,Mat &P,void *ctx) {

	// grab brinkman
	brinkman *par = static_cast<user_context*>(ctx)->prob;

	// get DM and system size
	DM sys,vel,pre;
	KSPGetDM(ksp,&sys);



	// get coarser brinkman
	brinkman *coarse = par->coarsen(sys);
	DMCompositeGetEntries(sys,&vel,&pre);

	DMCreateMatrix(pre,&P);

	// get element range
	int epai,epbi,epaj,epbj;
	brinkman::dm_fe_range(pre,par->po-1,epai,epbi,epaj,epbj);

	for (int ej=epaj; ej<epbj; ej++) for (int ei=epai; ei<epbi; ei++) {

		element ep(coarse->gr_p(ei,ej));
		coarse->add_P_mass(ep,pre,P);
	}

	MatAssemblyBegin(P,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(P,MAT_FINAL_ASSEMBLY);

	return 0;
}

/**
 * Calculate the stiffness matrix of the pressure system
 */
PetscErrorCode petsc_solve::calculate_p_stiff(KSP ksp,Mat &P,void *ctx) {

	// grab brinkman
	brinkman *par = static_cast<user_context*>(ctx)->prob;

	// get DM and system size
	DM sys,vel,pre;
	KSPGetDM(ksp,&sys);

	// get coarser brinkman
	brinkman *coarse = par->coarsen(sys);
	DMCompositeGetEntries(sys,&vel,&pre);

	DMCreateMatrix(pre,&P);

	// get element range
	int epai,epbi,epaj,epbj;
	brinkman::dm_fe_range(pre,par->po-1,epai,epbi,epaj,epbj);

	for (int ej=epaj; ej<epbj; ej++) for (int ei=epai; ei<epbi; ei++) {

		element ep(coarse->gr_p(ei,ej));
		coarse->add_P_stiff(ep,pre,P);
	}

	MatAssemblyBegin(P,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(P,MAT_FINAL_ASSEMBLY);

	return 0;
}

/**
 * Calculate the velocity operator
 */
void petsc_solve::calculate_A_vel(DM da_vel,Mat A,brinkman *par) {

	// get coarser brinkman + degrees of freedom
	brinkman *coarse = par->coarsen(da_vel);
	int &lu(coarse->lu),&lp(coarse->lp);

	// preallocation
	MatMPIAIJSetPreallocation(A,8*lu,NULL,8*lu,NULL);
	MatSeqAIJSetPreallocation(A,8*lu,NULL);
	MatSetUp(A);

	// get element range
	int euai,eubi,euaj,eubj;
	brinkman::dm_fe_range(da_vel,par->po,  euai,eubi,euaj,eubj);

	// loop through element grid
	for (int ej=euaj; ej<eubj; ej++) for (int ei=euai; ei<eubi; ei++) {

		// grab elements
		element eu(coarse->gr_u(ei,ej));

		// dump in velocities
		coarse->add_Auu(eu,da_vel,A);
		coarse->add_Avv(eu,da_vel,A);
	}

	// assemble
	MatAssemblyBegin(A,MAT_FLUSH_ASSEMBLY);
	MatAssemblyEnd(A,MAT_FLUSH_ASSEMBLY);

	// loop through element grid
	for (int ej=euaj; ej<eubj; ej++) for (int ei=euai; ei<eubi; ei++) {

		// grab elements
		element eu(coarse->gr_u(ei,ej));

		// place ones on velocity diagonal
		coarse->fix_AUU(eu,da_vel,A);
	}

	// finish up with assembly
	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
}

/**
 * Calculate the system operator
 */
void petsc_solve::calculate_A_sys(DM da_sys,Mat A,brinkman *par) {

	// get coarser brinkman + degrees of freedom
	brinkman *coarse = par->coarsen(da_sys);
	int &lu(coarse->lu),&lp(coarse->lp);

	// preallocation
	MatMPIAIJSetPreallocation(A,4*(2*lu+lp),NULL,4*(2*lu+lp),NULL);
	MatSeqAIJSetPreallocation(A,4*(2*lu+lp),NULL);
	MatSetUp(A);

	// get individual DMDAs
	DM da_vel,da_pre;
	DMCompositeGetEntries(da_sys,&da_vel,&da_pre);

	// get element range
	int euai,eubi,euaj,eubj;
	int epai,epbi,epaj,epbj;
	brinkman::dm_fe_range(da_vel,par->po,  euai,eubi,euaj,eubj);
	brinkman::dm_fe_range(da_pre,par->po-1,epai,epbi,epaj,epbj);

	// get index sets
	IS *lis,*gis;
	DMCompositeGetLocalISs( da_sys,&lis);
	DMCompositeGetGlobalISs(da_sys,&gis);

	// get submats
	Mat Auu,Aup,Apu,App;
	MatGetLocalSubMatrix(A,lis[0],lis[0],&Auu);
	MatGetLocalSubMatrix(A,lis[0],lis[1],&Aup);
	MatGetLocalSubMatrix(A,lis[1],lis[0],&Apu);
	MatGetLocalSubMatrix(A,lis[1],lis[1],&App);
	
	// loop through element grid
	for (int ej=euaj; ej<eubj; ej++) for (int ei=euai; ei<eubi; ei++) {

		// grab elements
		element eu(coarse->gr_u(ei,ej)), ep(coarse->gr_p(ei,ej));

		// dump in velocities
		coarse->add_Auu(eu,da_vel,Auu);
		coarse->add_Avv(eu,da_vel,Auu);

		// dump in vel-pre cross terms
		coarse->add_Auvp(eu,ep,da_vel,da_pre,Aup,Apu);
	}

	// put submats back
	MatRestoreLocalSubMatrix(A,lis[0],lis[0],&Auu);
	MatRestoreLocalSubMatrix(A,lis[0],lis[1],&Aup);
	MatRestoreLocalSubMatrix(A,lis[1],lis[0],&Apu);
	MatRestoreLocalSubMatrix(A,lis[1],lis[1],&App);

	// assemble
	MatAssemblyBegin(A,MAT_FLUSH_ASSEMBLY);
	MatAssemblyEnd(A,MAT_FLUSH_ASSEMBLY);

	// pull out diagonal submats
	MatGetLocalSubMatrix(A,lis[0],lis[0],&Auu);
	MatGetLocalSubMatrix(A,lis[1],lis[1],&App);

	// loop through element grid
	for (int ej=euaj; ej<eubj; ej++) for (int ei=euai; ei<eubi; ei++) {

		// grab elements
		element eu(coarse->gr_u(ei,ej)), ep(coarse->gr_p(ei,ej));

		// grab local pressure indicies
		std::vector<int> ip(coarse->lp);
		coarse->get_p_local_indices(ep,da_pre,ip.data());

		// place zeros on pressure diagonal
		for (int i:ip) MatSetValueLocal(App,i,i,0,INSERT_VALUES);

		// place ones on velocity diagonal
		coarse->fix_AUU(eu,da_vel,Auu);
	}

	// put submats back
	MatRestoreLocalSubMatrix(A,lis[1],lis[1],&App);
	MatRestoreLocalSubMatrix(A,lis[0],lis[0],&Auu);

	// finish up with assembly
	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

	// clean up index sets
	PetscFree(gis);
	PetscFree(lis);
}

/**
 * Calculate the system operator
 */
PetscErrorCode petsc_solve::calculate_A(KSP ksp,Mat A,Mat P,void *ctx) {

	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	// grab brinkman
	brinkman *par = static_cast<user_context*>(ctx)->prob;

	// get DM and system size
	DM dm;
	KSPGetDM(ksp,&dm);

	switch (brinkman::dm_type(dm)) {
		case brinkman::dmgr::sys: calculate_A_sys(dm,P,par); break;
		case brinkman::dmgr::vel: calculate_A_vel(dm,P,par); break;
		case brinkman::dmgr::pre:
			bail(-1,"this shouldn't happen\n");
			break;
	}

	return 0;
}

/**
 * initialization function: sets all variables to zero
 */
void brinkman::reset() {
	for (int i=0;i<gu;i++) {ux[i]=uy[i]=0; if(i < gp) p[i]=0;}
}

/**
 * Read a problem type from a saved binary file
 *
 * @param[in] fn filename
 * @return the problem type contained in the file
 */
ptype brinkman::read_type(const char *fn) {

	// get problem type
	ptype pt;
	FILE *fb = fopen(fn,"rb");
	fread(&pt,sizeof(ptype),1,fb);
	fclose(fb);

	return pt;
}


/**
 * allocate a brinkman object from a binary save file whose
 * brinkman subclass is not known a priori
 *
 * @param[in] fn filename
 * @param[out] ts pointer to problem_info
 */
brinkman* brinkman::load_binary(const char *fn) {

	// use subclass-specific allocation
	brinkman *info;
	switch (read_type(fn)) {

		case ptype::two_spheres:
			return two_spheres::load_binary(fn);

		case ptype::one_sphere:
			return one_sphere::load_binary(fn);

		case ptype::driven_cyl:
			return driven_cyl::load_binary(fn);

		default: bail(-1,"bad problem type\n"); return NULL;
	}
}

// solution accessors

/**
 * Get solution value at arbitrary point in the domain
 * @tparam which solution vector
 * @param[in] x horizontal coord
 * @param[in] y vertical coord
 * @return the solution
 */
template<dcomplex* brinkman::*arr>
dcomplex brinkman::interp_sol(double x,double y) {

	// reference to master element
	master_quad &mast = arr==&brinkman::p ? mast_p : mast_u;

	// reference to grid
	grid &gr = arr==&brinkman::p ? gr_p : gr_u;

	// number of nodes in an element
	int ln = arr==&brinkman::p ? lp : lu;
	dcomplex* val = new dcomplex[ln],retval=0;

	// requested solution array
	dcomplex* sol = this->*arr;

	// grab element index and local coords
	int ei,ej;
	double xi,eta;

	// grab element if we got one in the grid
	if (dm->inv_coords(x,y,ei,ej,xi,eta)) {

		// fill vals array
		for(node n:gr(ei,ej)) val[n.l]=sol[n.g];

		// evaluate result
		retval = mast.eval<dcomplex>(xi,eta,val);

	// otherwise, a nan
	} else retval = std::nan("out of grid")*(1+I);

	// cleanup and return
	delete[] val;
	return retval;
}

// explicit instantiation
template dcomplex brinkman::interp_sol<&brinkman::ux>(double,double);
template dcomplex brinkman::interp_sol<&brinkman::uy>(double,double);
template dcomplex brinkman::interp_sol<&brinkman::p>(double,double);

/**
 * Get solution & derivatives value at arbitrary point in the domain
 *
 * @tparam which solution vector
 * @param[in] x horizontal coord
 * @param[in] y vertical coord
 * @param[out] solution value
 * @param[out] solution x-derivative value
 * @param[out] solution y-derivative value
 */
template<dcomplex* brinkman::*arr>
void brinkman::interp_sol(double x,double y,
		dcomplex &v,dcomplex &v_x,dcomplex &v_y) {

	// reference to master element
	master_quad &mast = arr==&brinkman::p ? mast_p : mast_u;

	// reference to grid
	grid &gr = arr==&brinkman::p ? gr_p : gr_u;

	// number of nodes in an element
	int ln = arr==&brinkman::p ? lp : lu;
	dcomplex* val = new dcomplex[ln];

	// requested solution array
	dcomplex* sol = this->*arr;

	// grab element index and local coords
	int ei,ej;
	double xi,eta;

	// grab element if we got one in the grid
	if (dm->inv_coords(x,y,ei,ej,xi,eta)) {

		element e = gr(ei,ej);

		// grab coord trafo
		trafo_base *ct = dm->alloc_element_geometry(e);

		// fill vals array
		for(node n:e) val[n.l]=sol[n.g];

		// evaluate result
		mast.eval<dcomplex>(*ct,xi,eta,val,v,v_x,v_y);

		// free coord trafo
		dm->free_element_geometry(ct);

	} else {

		// otherwise set to nan
		v = v_x = v_y = std::nan("out of grid")*(1+I);
	}

	// cleanup
	delete[] val;
}

// explicit instantiation
template void brinkman::interp_sol<&brinkman::ux>(
		double,double,dcomplex&,dcomplex&,dcomplex&);
template void brinkman::interp_sol<&brinkman::uy>(
		double,double,dcomplex&,dcomplex&,dcomplex&);
template void brinkman::interp_sol<&brinkman::p>(
		double,double,dcomplex&,dcomplex&,dcomplex&);

/**
 * get the vorticity at (x,y)
 *
 * @param[in] x horizontal coord
 * @param[in] y vertical coord
 * @return the vorticity at the provided coordinate
 */
dcomplex brinkman::get_w(double x,double y) {
	dcomplex u,u_x,u_y,v,v_x,v_y;
	interp_sol<&brinkman::ux>(x,y,u,u_x,u_y);
	interp_sol<&brinkman::uy>(x,y,v,v_x,v_y);
	return v_x - u_y;
}

/**
 * get the divergence at (x,y)
 *
 * @param[in] x horizontal coord
 * @param[in] y vertical coord
 * @return the divergence at the provided coordinate
 */
dcomplex brinkman::get_d(double x,double y) {

	static const double tol=1e-6;

	dcomplex u,u_x,u_y,v,v_x,v_y;
	interp_sol<&brinkman::ux>(x,y,u,u_x,u_y);
	interp_sol<&brinkman::uy>(x,y,v,v_x,v_y);

	if (rev && fabs(y)<tol && fabs(v)>tol)
		bail(-1,"nonzero radial vel on axis\n");

	return rev ? (u_x+v_y+v/y) : (u_x+v_y);
}

/**
 * Calculate the complex advective term u.grad(u)
 * for a Brinkman solution at the provided location
 *
 * @param[in] x horizontal coord
 * @param[in] y vertical coord
 * @return the complex advective term u.grad(u)
 */
vcomplex brinkman::u_grad_u(double x,double y) {

	// grab solution, take conjugates
	dcomplex u,cu_x,cu_y,v,cv_x,cv_y;
	interp_sol<&brinkman::ux>(x,y,u,cu_x,cu_y);
	interp_sol<&brinkman::uy>(x,y,v,cv_x,cv_y);
	cu_x = std::conj(cu_x); cu_y = std::conj(cu_y);
	cv_x = std::conj(cv_x); cv_y = std::conj(cv_y);

	// advectiver term
	return vcomplex(u*cu_x + v*cu_y,u*cv_x + v*cv_y);
}

/**
 * Calculate the complex advective term u.grad(u)
 * for a Brinkman solution at the provided location,
 * calculating derivatives according to the chosen element
 *
 * @param[in] e element to use for derivatives
 * @param[in] n node at which to calculate
 * @return the complex advective term u.grad(u)
 */
vcomplex brinkman::u_grad_u(element &e,node &n) {

	// get element geometry
	trafo_base *ct = dm->alloc_element_geometry(e);

	// grab solution on element
	dcomplex *ux_sol = new dcomplex[lu];
	dcomplex *uy_sol = new dcomplex[lu];
	for(node nn:e) {ux_sol[nn.l]=ux[nn.g]; uy_sol[nn.l]=uy[nn.g];}

	// get point
	double xi,eta;
	mast_u.point(n.l,xi,eta);

	// grab solution, take conjugates
	dcomplex u,cu_x,cu_y,v,cv_x,cv_y;
	mast_u.eval(*ct,xi,eta,ux_sol,u,cu_x,cu_y);
	mast_u.eval(*ct,xi,eta,uy_sol,v,cv_x,cv_y);
	cu_x = std::conj(cu_x); cu_y = std::conj(cu_y);
	cv_x = std::conj(cv_x); cv_y = std::conj(cv_y);

	dm->free_element_geometry(ct);

	return vcomplex(u*cu_x + v*cu_y,u*cv_x + v*cv_y);
}

/**
 * Read out the petsc solution into the ux, uy and p arrays
 */
void brinkman::read(petsc_solve &solver,Vec X) {

	// by default, use solution from solver
	if (X==NULL) KSPGetSolution(solver.ksp,&X);

	// check how many processors we have
	int size;
	MPI_Comm_size(PETSC_COMM_WORLD,&size);

	// grab global indices
	IS *gis;
	DM sys;
	KSPGetDM(solver.ksp,&sys);
	DMCompositeGetGlobalISs(sys,&gis);

	// grab each DMDA
	DM vel,pre;
	DMCompositeGetEntries(sys,&vel,&pre);

	// grab separate velocity and pressure vectors
	Vec Xug,Xpg;
	if (size == 1) {

		// if we only have one proc, just grab solution vector

		// nest vs subvector?
		VecGetSubVector(X,gis[0],&Xug);
		VecGetSubVector(X,gis[1],&Xpg);

	} else {

		// otherwise, scatter to all procs
		
		// vel and pres parts
		Vec Xu,Xp,Xun,Xpn;

		// grab global vectors (these are ordered globally)
		VecGetSubVector(X,gis[0],&Xu);
		VecGetSubVector(X,gis[1],&Xp);

		// create vectors for natural ordering
		DMDACreateNaturalVector(vel,&Xun);
		DMDACreateNaturalVector(pre,&Xpn);


		// scatter into natural odering
		DMDAGlobalToNaturalBegin(vel,Xu,INSERT_VALUES,Xun);
		DMDAGlobalToNaturalBegin(pre,Xp,INSERT_VALUES,Xpn);
		DMDAGlobalToNaturalEnd(  pre,Xp,INSERT_VALUES,Xpn);
		DMDAGlobalToNaturalEnd(  vel,Xu,INSERT_VALUES,Xun);

		// now set up scatters to all procs
		VecScatter scat_u,scat_p;
		VecScatterCreateToAll(Xun,&scat_u,&Xug);
		VecScatterCreateToAll(Xpn,&scat_p,&Xpg);

		VecScatterBegin(scat_u,Xun,Xug,INSERT_VALUES,SCATTER_FORWARD);
		VecScatterBegin(scat_p,Xpn,Xpg,INSERT_VALUES,SCATTER_FORWARD);
		VecScatterEnd(  scat_p,Xpn,Xpg,INSERT_VALUES,SCATTER_FORWARD);
		VecScatterEnd(  scat_u,Xun,Xug,INSERT_VALUES,SCATTER_FORWARD);


		// clean up
		VecDestroy(&Xun);
		VecDestroy(&Xpn);
		VecScatterDestroy(&scat_u);
		VecScatterDestroy(&scat_p);
		//DMCompositeRestoreAccess(da_sys,X,&Xu,&Xp);
		VecRestoreSubVector(X,gis[0],&Xu);
		VecRestoreSubVector(X,gis[1],&Xp);
	}

	// array values
	dcomplex *xu,*xp;
	VecGetArray(Xug,&xu);
	VecGetArray(Xpg,&xp);

	// copy values into ux/uy/p
	for (int i=0;i<gu;i++) {
		ux[i] = xu[2*i];
		uy[i] = xu[2*i+1];
	}

	for (int i=0;i<gp;i++) p[i] = xp[i];

	VecRestoreArray(Xug,&xu);
	VecRestoreArray(Xpg,&xp);

	if (size == 1) {
		VecRestoreSubVector(X,gis[0],&Xug);	
		VecRestoreSubVector(X,gis[1],&Xpg);	
	} else {

		VecDestroy(&Xug);
		VecDestroy(&Xpg);
	}

	ISDestroy(&gis[0]);
	ISDestroy(&gis[1]);
	PetscFree(gis);
}

/**
 * generalized gp quad output
 *
 * @tparam fld which field to output (0-pressure, 1-x velocity,
 *   2-y velocity)
 * @tparam real whether to output the real part of the field
 * @param[in] fn string containing filename to write to
 */
template<int fld,bool real> void brinkman::gp_quad(const char* fn) {

	if (!gru()) return;

	// whether this is the pressure
	static const bool pp=fld==0;
	// whether this is the y-velocity
	static const bool  y=fld==1;

	FILE *fh = fopen(fn,"w");

	// loop through pressure or vel grid as needed
	for (element e: pp?gr_p:gr_u) {

		// grab geometry and allocate data
		trafo_base *ct = dm->alloc_element_geometry(e);
		double *dat = new double[e.n_nds];

		// read in value from field
		for (node n:e) {
			dcomplex arg = pp?p[n.g]:(y?uy[n.g]:ux[n.g]); 
			dat[n.l] = real?std::real(arg):std::imag(arg);
		}

		// output data and deallocate
		(pp?mast_p:mast_u).gp_fill(*ct,dat,fh);
		delete[] dat;
		dm->free_element_geometry(ct);
	}

	fclose(fh);
}

/**
 * generalized gp lines output
 *
 * @tparam fld which field to output (0-pressure, 1-x velocity,
 *   2-y velocity)
 * @tparam real whether to output the real part of the field
 * @param[in] fn string containing filename to write to
 */
template<int fld,bool real> void brinkman::gp_line(const char *fn) {

	if (!gru()) return;

	// whether this is the pressure
	static const bool pp=fld==0;
	// whether this is the y-velocity
	static const bool  y=fld==1;

	FILE *fh = fopen(fn,"w");

	// loop through pressure or vel grid as needed
	for(element e : pp?gr_p:gr_u) {

		// grab geometry and allocate data
		trafo_base *ct = dm->alloc_element_geometry(e);
		double *dat = new double[e.n_nds];

		// read in values from field
		for (node n:e) {
			dcomplex arg = pp?p[n.g]:(y?uy[n.g]:ux[n.g]);
			dat[n.l] = real?std::real(arg):std::imag(arg);
		}

		// output data and deallocate
		(pp?mast_p:mast_u).gp_lines(*ct,dat,fh);
		delete[] dat;
		dm->free_element_geometry(ct);
	}

	fclose(fh);
}

// explicit instantiation
template void brinkman::gp_quad<0,true>(const char*);
template void brinkman::gp_quad<1,true>(const char*);
template void brinkman::gp_quad<2,true>(const char*);
template void brinkman::gp_quad<0,false>(const char*);
template void brinkman::gp_quad<1,false>(const char*);
template void brinkman::gp_quad<2,false>(const char*);
template void brinkman::gp_line<0,true>(const char*);
template void brinkman::gp_line<1,true>(const char*);
template void brinkman::gp_line<2,true>(const char*);
template void brinkman::gp_line<0,false>(const char*);
template void brinkman::gp_line<1,false>(const char*);
template void brinkman::gp_line<2,false>(const char*);
