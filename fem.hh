///////////////////////////////////////////////////////////////////////////////
// fem.hh
//
// a collection of structs for implementing the finite element method over
// logically rectangular grids on arbitrary geometries via coordinate
// transformation functions
///////////////////////////////////////////////////////////////////////////////

#ifndef FEM_HH
#define FEM_HH

#include <algorithm>
#include <vector>
#include <cmath>
#include <cstdio>
#include "petscksp.h"

/**
 * represents a single global node of a finite element
 */
struct gnode {
	/** local and global indices referring to this node */
	int g;
	/** global col/row */
	int gj,gi;

	/** empty constructor sets all to zero */
	gnode():g(0),gj(0),gi(0){}
	/**
	 * constructor setting l
	 *
	 * @param[in] l local node index
	 * @param[in] li local horizontal index
	 * @param[in] lj local vertical index
	 */
	gnode(int g_,int gi_,int gj_):g(g_),gj(gj_),gi(gi_){}
};

/**
 * represents a single node of a finite element
 */
struct node {
	/** local and global indices referring to this node */
	int l,g;
	/** local col/row */
	int lj,li;
	/** global col/row */
	int gj,gi;

	/** empty constructor sets all to zero */
	node():l(0),g(0),lj(0),li(0),gj(0),gi(0){}
	/**
	 * constructor setting l
	 *
	 * @param[in] l local node index
	 * @param[in] li local horizontal index
	 * @param[in] lj local vertical index
	 */
	node(int l_,int li_,int lj_):l(l_),g(0),lj(lj_),li(li_),gj(0),gi(0){}
	node(int l_,int li_,int lj_,int g_,int gi_,int gj_):
		l(l_),g(g_),lj(lj_),li(li_),gj(gj_),gi(gi_) {}

};

/**
 * use node pairs to find neighbors, e.g. for
 * inverting local stiffness mat for additive schwarz
 */
typedef std::pair<node,node> node_pair;

/**
 * represents a quad finite element in a logical grid
 *
 * an element of order p contains (p+1)^2 nodes, with
 * local indices as e.g. for p=3,
 *
 * 12 13 14 15
 *  8  9 10 11
 *  4  5  6  7
 *  0  1  2  3
 *
 */
struct element {

	// iterator

	/**
	 * iterator for looping through nodes in an element,
	 * inheriting from std::input_iterator_tag
	 */
	struct iterator: public std::iterator<
		std::input_iterator_tag,node,int,const node*,node> {

		/** polynomial order */
		int p;
		/** actual nodes in grid dirs (ignoring periodicity) */
		int am,an;
		/** grid periodicity in horiz, vert direction */
		bool hp,vp;
		/** corner indices */
		int ai,aj;
		/** node */
		node n;

		/** unique nodes in horiz direction (including periodicity) */
		int um() {return hp?(am-1):am;}
		/** unique nodes in vert direction (including periodicity) */
		int un() {return vp?(an-1):an;}

		void check_index() {

			// allow "end" node
			if (n.li==0 && n.lj==p+1) return;

			// check horizontal index
			if (n.gi<0 || n.gi>=am) {
				fprintf(stderr,
						"global horiz index %d out of range [0,%d) \n",n.gi,am);
				exit(-1);
			}
			if (n.gi == um()) {
				n.gi=0;
			}

			// check vertical index
			if (n.gj<0 || n.gj>=an) {
				fprintf(stderr,
						"global vert index %d out of range [0,%d) \n",n.gj,an);
				exit(-1);
			}
			if (n.gj == un()) {
				n.gj=0;
			}
		}

		/** set global index based on local vals */
		void set_global();

		iterator(int p_,int me,int ne,bool h_prd,bool v_prd,
				int ai_,int aj_,int nl,int ni,int nj):
			p(p_),am(me*p+1),an(ne*p+1),hp(h_prd),vp(v_prd),
			ai(ai_),aj(aj_),n(nl,ni,nj) {}

		iterator(int p,int me,int ne,bool h_prd,bool v_prd,int ai,int aj):
			iterator(p,me,ne,h_prd,v_prd,ai,aj,0,0,0) {}

		iterator(element &e,int nl,int ni,int nj):
			iterator(e.p,e.m,e.n,e.h_prd,e.v_prd,e.i*e.p,e.j*e.p,nl,ni,nj) {
				set_global();
		}

		iterator(element &e): iterator(e,0,0,0) {set_global();}


		/** pre-increment */
		iterator& operator++();
		/** post-increment */
		iterator operator++(int) {iterator ret=*this; ++(*this); return ret;}
		/** equality is based on local nodes */
		bool operator==(iterator other) const {return n.l==other.n.l;}
		/** inequality */
		bool operator!=(iterator other) const {return !(*this==other);}
		/** return the node when dereferenced */
		reference operator*() const {return n;}
	};

	/** returns an iterator set at the first node */
	iterator begin() {
		return in_grid()?iterator(*this):end();
	}
	/** returns an iterator set one past the final node */
	iterator end() {return iterator(*this,(p+1)*(p+1),0,p+1);}

	// element data members

	/** polynomial order */
	int p;
	/** i,j coords of element in grid */
	int i,j;
	/** grid extent, horizontally and vertically */
	int m,n;
	/** grid periodicity, horizontal and vertical */
	bool h_prd,v_prd;
	/** number of nodes */
	int n_nds;
	/** linear 1D index of element in grid */
	int ind() const {return i+j*m;}

	/** empty constructor */
	element(int p_,int i_,int j_,int m_,int n_,bool h_prd_,bool v_prd_):
		p(p_),i(i_),j(j_),m(m_),n(n_),h_prd(h_prd_),v_prd(v_prd_),
		n_nds((p+1)*(p+1)) {}

	element():element(0,0,0,0,0,false,false){}

	/**
	 * constructs an element of specified index
	 * in a grid with the specified limits
	 *
	 * @param[in] p polynomial order
	 * @param[in] i horizontal element index
	 * @param[in] j vertical element index
	 * @param[in] m number of elements in horiz direction
	 * @param[in] n number of elements in vertical direction
	 */
	element(int p,int i,int j,int m,int n):
		element(p,i,j,m,n,false,false) {}



	/** elements are equal if they have the same row/col index */
	bool operator==(element other) const {return i==other.i&&j==other.j;}
	/** elements are not equal if they have a different index */
	bool operator!=(element other) const {return !(*this==other);}
	/** less than is based on global row-major ordering */
	bool operator<(element other) const {return j<other.j || (i<other.i&&j==other.j);}
	/** less than / equals */
	bool operator<=(element other) const {return (*this)==other||(*this)<other;}
	/** greater than */
	bool operator>(element other) const {return !((*this)<=other);}
	/** greater than / equals */
	bool operator>=(element other) const {return !((*this)<other);}

	/**
	 * returns an element located a provided number of
	 * horizontal and vertical indices away
	 *
	 * @param[in] di difference in horizontal index
	 * @param[in] dj difference in vertical index
	 */
	element neigh(int di,int dj) {return element(p,i+di,j+dj,m,n,h_prd,v_prd);}

	/** returns true if this element is in the grid */
	bool in_grid() {return i>=0&&i<m&&j>=0&&j<n;}

	/** get vector of all neighboring elements in grid */
	void get_neighbors(std::vector<element> &es);

	/** get list of nodes shared with another element */
	void get_shared_nodes(element other,std::vector<node_pair> &shared);
};

/**
 * represents a cartesian grid of quad elements
 */
struct grid {

	/**
	 * iterator implementing std::random_access_iterator_tag,
	 * allowing for the use of OpenMP parallelization given
	 * a range-based for loop
	 */
	struct iterator: public std::iterator<
		std::random_access_iterator_tag,element,int,const element*,element> {

		/** pointer to parent grid */
		grid *g;
		/** current element */
		element e;

		/** empty constructor */
		iterator(): g(NULL) {}
		/** constructor which inializes iterator to first element */
		iterator(grid &g_) : g(&g_), e(g->p,0,0,g->m,g->n,g->h_prd,g->v_prd) {}
		/** constructor which inializes iterator to specified element */
		iterator(grid &g_,int ie,int je):
			g(&g_),e(g->p,ie,je,g->m,g->n,g->h_prd,g->v_prd) {}

		/** pre-increment */
		iterator& operator++();
		/** post-increment */
		iterator operator++(int) {iterator ret=*this;++(*this);return ret;}
		/** pre-decrement */
		iterator& operator--();
		/** post-decrement */
		iterator operator--(int) {iterator ret=*this; --(*this); return ret;}

		// all comparisons are based on elements

		/** equals */
		bool operator==(iterator other) const {return e==other.e;}
		/** not equals */
		bool operator!=(iterator other) const {return e!=other.e;}
		/** less than */
		bool operator<(iterator other) const {return e<other.e;}
		/** greater than */
		bool operator>(iterator other) const {return e>other.e;}
		/** less than / equals */
		bool operator<=(iterator other) const {return e<=other.e;}
		/** greater than / equals */
		bool operator>=(iterator other) const {return e>=other.e;}

		/** addition */
		iterator operator+(int n) const;
		/** subtraction */
		iterator operator-(int n) const {return (*this)+(-n);}
		/** iterator difference (global index) */
		int operator-(const iterator other) const {return e.ind() - other.e.ind();}
		/** add/equals */
		iterator operator+=(int n);
		/** subtraction/equals */
		iterator operator-=(int n) {return (*this) += (-n);}
		/** return element when dereferenced */
		reference operator*() const {return e;}

	};

	/** polynomial order */
	int p;
	/** horizontal/vertical # of elements */
	int m,n;
	/** periodicity in horiz / vert directions */
	bool h_prd,v_prd;
	/** number of elements */
	int n_els;
	/** number of points in horiz/vertical directions */
	int mp,np;
	/** number of dof */
	int n_dof;

	grid() {}
	/**
	 * Constructor which assembles a finite element of specified
	 * order and number of elements
	 *
	 * @param[in] p polynomial order
	 * @param[in] m number of horizontal elements
	 * @param[in] n number of vertical elements
	 * @param[in] h_prd horizontal periodicity
	 * @param[in] v_prd vertical periodicity
	 */
	grid(int p_,int m_,int n_,bool hp,bool vp):
		p(p_),m(m_),n(n_),h_prd(hp),v_prd(vp),n_els(m*n),
		mp(m*p+(hp?0:1)),np(n*p+(vp?0:1)),n_dof(mp*np) {}
	/**
	 * Constructor which assembles a finite element of specified
	 * order and number of elements
	 *
	 * @param[in] p polynomial order
	 * @param[in] m number of horizontal elements
	 * @param[in] n number of vertical elements
	 */
	grid(int p,int m,int n):grid(p,m,n,false,false) {}

	/**
	 * Parentheses operator allows access to an element at
	 * the provided index pair
	 *
	 * @param[in] i horizontal index
	 * @param[in] j vertical index
	 * @return the element at (i,j) in the grid
	 */
	element operator()(int i,int j){return element(p,i,j,m,n,h_prd,v_prd);}

	/** iterator pointing to first element in grid */
	iterator begin() {return iterator(*this);}
	/** iterator pointing to one after the last element in grid */
	iterator end() {return iterator(*this,0,n);}

	void mg_grid_info(int level,grid &fine,grid &coar) {

		int ffac = 1<<level;
		int cfac = ffac<<1;

		fine = grid(p,m/ffac,n/ffac);
		coar = grid(p,m/cfac,n/cfac);
	}

	/**
	 */
	void fill_petsc_interp(int level,double *interp,
			Mat &I,int row_off=0,int col_off=0) {

		grid fine,coar;
		mg_grid_info(level,fine,coar);


		// multiplicity
		int *mult = new int[fine.n_dof];
		for(int i=0;i<fine.n_dof;i++) mult[i]=0;

		// look at the corresponding 2x2 grid
		// of fine nodes corresponding to this coarse node
		grid lf(coar.p,2,2);

		// nodes per element
		int nds = coar(0,0).n_nds;

		// for each element, store the (p+1)^2 x (2p+1)^2 interpolation matrix
		// ...or WOULD do that, but each is just the passed in interp matrix

		// row/col inds
		int *row_inds = new int[lf.n_dof];
		int *col_inds = new int[nds];

		// marker array for not double counting
		bool *visit = new bool[lf.n_dof];

		// go through each coarse element to count up multiplicites
		for(element ce:coar) {


			// only visit each fine grid over this patch once
			for(int i=0;i<lf.n_dof;i++) visit[i]=false;

			// go through each local fine element
			for(element lfe:lf) {

				// get equivalent global fine element
				element gfe = fine(2*ce.i + lfe.i,2*ce.j + lfe.j);

				// grab global iterator
				element::iterator g_it = gfe.begin();

				// loop over local fine nodes
				for (node fn: lfe) {

					// loop through coarse nodes if haven't
					// been here before
					if (!visit[fn.g]) mult[(*g_it).g]++;

					// mark that we've visited
					visit[fn.g]=true;

					// increment iterator
					g_it++;
				}
			}
		}

		// correct for multiplicity
		std::complex<double> *vals = new std::complex<double>[nds * lf.n_dof];

		// go through each coarse element
		for(element ce:coar) {

			// set col inds
			for (node n: ce) col_inds[n.l] = n.g + col_off;

			// copy in master interp
			for(int i=0;i<nds*lf.n_dof;i++) vals[i] = interp[i];

			// visit each fine node once
			for(int i=0;i<lf.n_dof;i++) visit[i]=false;

			// go through each local fine element
			for(element lfe:lf) {

				// get equivalent global fine element
				element gfe = fine(2*ce.i + lfe.i,2*ce.j + lfe.j);

				// grab global iterator
				element::iterator g_it = gfe.begin();

				// loop over local fine nodes
				for (node fn: lfe) {

					// if we haven't been here before...
					if (!visit[fn.g]) {

						// global fine node from iterator
						int gf = (*g_it).g;

						// store row index
						row_inds[fn.g] = gf + row_off;

						// get multiplicity
						int mul = mult[gf];

						// correct for each coarse node
						for (node cn:ce) {

							// global coarse node
							int gc = cn.g;

							// alter value
							vals[cn.l + fn.g*nds] /= mul;

						}

						// mark that we've visited
						visit[fn.g]=true;
					}

					// increment iterator
					g_it++;
				}
			}

			// add this to itnerp matrix
			MatSetValues(I,lf.n_dof,row_inds,nds,col_inds,vals,ADD_VALUES);
		}

		delete[] mult;
		delete[] visit;
		delete[] row_inds;
		delete[] col_inds;
	}

	void free_interp(double *(&interp)) {delete[] interp; interp=NULL;}
};

/**
 * represents a cartesian grid of global nodes
 */
struct gnode_grid {

	/**
	 * iterator implementing std::random_access_iterator_tag,
	 * allowing for the use of OpenMP parallelization given
	 * a range-based for loop
	 */
	struct iterator: public std::iterator<
		std::random_access_iterator_tag,gnode,int,const gnode*,gnode> {

		/** pointer to parent grid */
		gnode_grid *g;
		/** current node */
		gnode n;

		/** empty constructor */
		iterator(): g(NULL) {}
		/** constructor which inializes iterator to first element */
		iterator(gnode_grid &g_) : g(&g_) {}
		/** constructor which inializes iterator to specified element */
		iterator(gnode_grid &g_,int i,int j):g(&g_),n(i+j*g->m,i,j) {}

		/** pre-increment */
		iterator& operator++();
		/** post-increment */
		iterator operator++(int) {iterator ret=*this;++(*this);return ret;}
		/** pre-decrement */
		iterator& operator--();
		/** post-decrement */
		iterator operator--(int) {iterator ret=*this; --(*this); return ret;}

		// all comparisons are based on elements

		/** equals */
		bool operator==(iterator other) const {return n.g==other.n.g;}
		/** not equals */
		bool operator!=(iterator other) const {return n.g!=other.n.g;}
		/** less than */
		bool operator<(iterator other) const {return n.g<other.n.g;}
		/** greater than */
		bool operator>(iterator other) const {return n.g>other.n.g;}
		/** less than / equals */
		bool operator<=(iterator other) const {return n.g<=other.n.g;}
		/** greater than / equals */
		bool operator>=(iterator other) const {return n.g>=other.n.g;}

		/** addition */
		iterator operator+(int d) const;
		/** subtraction */
		iterator operator-(int d) const {return (*this)+(-d);}
		/** iterator difference (global index) */
		int operator-(const iterator other) const {return n.g - other.n.g;}
		/** add/equals */
		iterator operator+=(int d);
		/** subtraction/equals */
		iterator operator-=(int d) {return (*this) += (-d);}
		/** return element when dereferenced */
		reference operator*() const {return n;}

	};

	/** horizontal/vertical # of nodes */
	int m,n;
	/** number of elements */
	int n_nds;

	gnode_grid() {}
	/**
	 * Constructor which assembles a grid of nodes
	 *
	 * @param[in] m number of horizontal nodes
	 * @param[in] n number of vertical nodes
	 */
	gnode_grid(int m_,int n_):m(m_),n(n_),n_nds(m*n) {}


	/** iterator pointing to first element in grid */
	iterator begin() {return iterator(*this);}
	/** iterator pointing to one after the last element in grid */
	iterator end() {return iterator(*this,0,n);}

};

/** non-member addition */
grid::iterator operator+(int n,const grid::iterator it);
/** non-member subtraction */
grid::iterator operator-(int n,const grid::iterator it);
/** non-member addition */
gnode_grid::iterator operator+(int n,const gnode_grid::iterator it);
/** non-member subtraction */
gnode_grid::iterator operator-(int n,const gnode_grid::iterator it);

// coordinate transforms and domains

/**
 * abstract class for an element-specific coordinate transformation
 */
struct trafo_base {
	/** get real coordinates (x,y) at master coordinates (xi,eta) */
	virtual void xi2x(double xi,double et,double &x,double &y) const=0;
	/** get determinant of the Jacobean at master coords (xi,eta) */
	virtual double Jdet(double xi,double et) const=0;
	/** get inverse Jacobean at (xi,eta) */
	virtual void Jinv(double xi,double et,
			double &xi_x,double &et_x,double &xi_y,double &et_y) const=0;
	virtual void Jinv2(double xi,double et,
		double &xi_xx, double &xi_xy, double &xi_yy,
		double &et_xx, double &et_xy, double &et_yy) const=0;
	/** get inverse Jacobean derivatives at (xi,eta) */
	virtual void Jinv_det(double xi,double et,
			double &xi_x_xi,double &et_x_xi,double &xi_y_xi,double &et_y_xi,
			double &xi_x_et,double &et_x_et,double &xi_y_et,double &et_y_et) const{
		fputs("trafo_base::Jinv_det not implemented by default\n",stderr);
		exit(-1);
	}
	/** get coordinates (xi,eta) at coordinate (x,y) */
	virtual void x2xi(double x,double y,double &xi,double &eta) const {
		fputs("trafo_base::x2xi not implemented by default\n",stderr);
		exit(-1);
	}
};


/**
 * abstract interface for domain -> element-wise coordinate
 * tranformation information
 */
struct domain {

	/**
	 * returns a reference to a newly allocated trafo_base
	 * object corresponding to the provided element
	 *
	 * @param[in] e the element to pull a coordinate transform for
	 * @return a pointer to a new trafo_base object
	 */
	virtual trafo_base* alloc_element_geometry(element &e)=0;

	/**
	 * deallocates the trafo_base pointed to by the provided pointer
	 *
	 * @param[in] ct pointer to a trafo_base object
	 * @param[out] ct the same pointer, now set to NULL
	 */
	void free_element_geometry(trafo_base *(&ct))
		{if (ct!=NULL) delete ct; ct=NULL;}

	/**
	 * returns the element indices and intra-element coordinates
	 * corresponding to the provided top level global coordinate
	 *
	 * @param[in] x horiz global coord
	 * @param[in] y vert global coord
	 * @param[out] ei horiz element index
	 * @param[out] ej vert element index
	 * @param[out] xi in-element horiz coord
	 * @param[out] eta in-element vert coord
	 */
	virtual bool inv_coords(double x,double y,
		int &ei,int &ej,double &xi,double &eta) const=0;
};


/**
 * struct which relates an element to a Cartesian grid
 *
 * optionally, represents two grids glued together 
 * at an intermediate x-value
 */
struct cartesian_grid: public domain {

	/** horizontal and vertical grid starting point */
	double ax,ay;
	/** horizontal and vertical grid ending point */
	double bx,by;
	/** width/height of an element */
	double ex,ey;

	/**
	 * constructs a grid of the specified number of elements with
	 * the provided extent in cartesian space
	 *
	 * @param m horizontal number of elements
	 * @param n vertical number of elements
	 * @param ax left-edge cartesian coordinate
	 * @param ay bottom-edge cartesian coordinate
	 * @param bx right-edge cartesian coordinate
	 * @param by top-edge cartesian coordinate
	 */
	cartesian_grid(int m,int n,double ax_,double bx_,double ay_,double by_):
		ax(ax_),ay(ay_),bx(bx_),by(by_),ex((bx_-ax_)/m),ey((by_-ay_)/n) {
	}

	/**
	 * coordinate transormation to a square piece of a cartesian domain
	 */
	struct linear_geom: public trafo_base {

		/** element from which coords are transformed */
		element &e;
		/** left/bottom cartesian coordinate */
		double aex,aey;
		/** width/height of this patch of domain */
		double dx,dy;

		/** constructor explicitly setting each class member */
		linear_geom(element e_,double aex_,double aey_,double dx_,double dy_):
			e(e_),aex(aex_),aey(aey_),dx(dx_),dy(dy_) {
		}

		/** constructor pulling information from parent cartesian_grid */
		linear_geom(cartesian_grid &p,element &e_):
			e(e_),aex(p.ax+e.i*p.ex),aey(p.ay+e.j*p.ey),dx(p.ex),dy(p.ey) {
			}

		/** get real coordinates (x,y) at master coordinates (xi,eta) */
		void xi2x(double xi,double et,double &x,double &y) const
			{x=aex+xi*dx;y=aey+et*dy;}
		/** get determinant of the Jacobean at master coords (xi,eta) */
		double Jdet(double xi,double et) const {return dx*dy;}
		/** get inverse Jacobean at (xi,eta) */
		void Jinv(double xi,double et,double &xi_x,
			double &et_x,double &xi_y,double &et_y) const
			{
				xi_x=1./dx;et_x=0;xi_y=0;et_y=1./dy;
			}
		void Jinv2(double xi,double et,
			double &xi_xx, double &xi_xy, double &xi_yy,
			double &et_xx, double &et_xy, double &et_yy) const {
			xi_xx=xi_xy=xi_yy=et_xx=et_xy=et_yy = 0;
		}
		void Jinv_det(double xi,double et,
			double &xi_x_xi,double &et_x_xi,double &xi_y_xi,double &et_y_xi,
			double &xi_x_et,double &et_x_et,double &xi_y_et,double &et_y_et) const{
			xi_x_xi=et_x_xi=xi_y_xi=et_y_xi=xi_x_et=et_x_et=xi_y_et=et_y_et=0;
		}
	};

	/**
	 * coordinate transform using this cartesian grid as an intermediate
	 * to another more complicated geometry
	 */
	struct trafo_geom : public trafo_base {

		/** linear geometry trafo */
		linear_geom p;
		/**
		 * another trafo detailing the relationship betweeen
		 * this 2D cartesian and another geometry
		 */
		trafo_base &ct;

		/** 
		 * constructor creating member linear_geom from 
		 * reference to parent cartesian_grid
		 */
		trafo_geom(cartesian_grid &g,element &e,trafo_base &ct_):p(g,e),ct(ct_) {}

		/** 
		 * constructor creating member linear_geom from copy
		 */
		trafo_geom(linear_geom p_,trafo_base &ct_):p(p_),ct(ct_) {}

		/** get real coordinates (x,y) at master coordinates (xi,eta) */
		void xi2x(double xi,double et,double &x,double &y) const;
		/** get determinant of the Jacobean at master coords (xi,eta) */
		double Jdet(double xi,double et) const;
		/** get inverse Jacobean at (xi,eta) */
		void Jinv(double xi,double et,double &xi_x,double &et_x,
				double &xi_y,double &et_y) const;
		void Jinv2(double xi,double et,
			double &xi_xx, double &xi_xy, double &xi_yy,
			double &et_xx, double &et_xy, double &et_yy) const;
		void Jinv_det(double xi,double et,
			double &xi_x_xi,double &et_x_xi,double &xi_y_xi,double &et_y_xi,
			double &xi_x_et,double &et_x_et,double &xi_y_et,double &et_y_et) const;
	};

	/** returns the linear geom at element e */
	virtual linear_geom element_geometry(element &e)
		{return linear_geom(*this,e);}
	/** allocates the linear geom at element e */
	virtual trafo_base* alloc_element_geometry(element &e)
		{ //printf("here: [%g,%g] x [%g,%g] (%g,%g)\n",ax,bx,ay,by,ex,ey);
			return new linear_geom(*this,e);}

	/**
	 * creates a compout coord transformation at element e
	 * going from element space to this intermediate grid
	 * to a space defined by the provided trafo_base
	 */
	trafo_geom* coord_trafo(element &e,trafo_base &t) {
		return new trafo_geom(element_geometry(e),t);}

	/**
	 * Returns the element index and master
	 * coodinates in grid/element space corresponding
	 * to the provided cartesian coordinate, trafo'd
	 * thru and intermeidate coordinate transform
	 *
	 * @param[in] ct intermediate transform
	 * @param[in] x cartesian (outer) coordinate
	 * @param[in] y cartesian (outer) coordinate
	 * @param[out] ei horizontal element index of (x,y)
	 * @param[out] ej vertical element index of (x,y)
	 * @param[out] xi horiz master coordinate w/in elem
	 * @param[out] eta vert master coordinate w/in elem
	 * @return true iff the provided point is in the grid
	 */
	virtual bool inter_inv_coords(const trafo_base &ct,double x,double y,
			int &ei,int &ej,double &xi,double &eta) const {

		// get intermediate coords
		double x_i,y_i;
		ct.x2xi(x,y,x_i,y_i);

		// if x or y is out of domain, return false
		if (x_i<ax || x_i>bx || y_i<ay || y_i>by) return false;

		// relate local coords to dist past "left" edge
		xi  = (x_i-ax)/ex;
		eta = (y_i-ay)/ey;

		// indices are integer part, so store and subtract
		ei = xi>0  ? (int) xi  : (((int) xi)  - 1);
		xi-=ei;

		ej = eta>0 ? (int) eta : (((int) eta) - 1);
		eta-=ej;

		// the point was in the domain
		return true;
	}
	
	/**
	 * Returns element index and local / master coordinates
	 * at a given cartesian coordinate
	 */
	virtual bool inv_coords(double x,double y,
			int &ei,int &ej,double &xi,double &eta) const {

		// if x or y is out of domain, return false
		if (x<ax || x>bx || y<ay || y>by) return false;

		// relate local coords to dist past "left" edge
		xi  = (x-ax)/ex;
		eta = (y-ay)/ey;

		// indices are integer part, so store and subtract
		ei = xi>0  ? (int) xi  : (((int) xi)  - 1);
		xi-=ei;

		ej = eta>0 ? (int) eta : (((int) eta) - 1);
		eta-=ej;

		// the point was in the domain
		return true;
	}
};

/**
 * structure representing two regular cartesian grids glued
 * together at an intermediate x-value
 */
struct xsplit_grid: public cartesian_grid {

	/** index of the grid split */
	int mi;
	/** "right"-most horizontal cartesian coord */
	double cx;
	/** element width in the second half of the split grid */
	double ex2;

	/**
	 * constructor
	 *
	 * @param[in] m1 number of columns of elements in the first (x) grid
	 * @param[in] m2 number of columns of elements in the second (x) grid
	 * @param[in] n  number of rows of elements in the (y) grid
	 * @param[in] ax left edge coordinate
	 * @param[in] bx intermediate coordinate
	 * @param[in] cx right edge coordinate
	 * @param[in] ay bottom edge coordinate
	 * @param[in] by top edge coordinate
	 */
	xsplit_grid(int m1,int m2,int n,double ax,double bx,double cx_,double ay,double by):
		cartesian_grid(m1,n,ax,bx,ay,by),mi(m1),cx(cx_),ex2((cx-bx)/m2) {}

	/** helper function returning the left x-coord of element i */
	double calc_ax(int i) {return ax + (i<mi?(i*ex):(mi*ex+(i-mi)*ex2));}
	/** helper function returning the width of element i */
	double calc_ex(int i) {return i<mi?ex:ex2;}

	/** returns the linear geom for element e */
	linear_geom element_geometry(element &e) {
		return linear_geom(e,calc_ax(e.i),ay+e.j*ey,calc_ex(e.i),ey);}
	/** allocates the linear geom for element e */
	trafo_base* alloc_element_geometry(element &e) {
		return new linear_geom(e,calc_ax(e.i),ay+e.j*ey,calc_ex(e.i),ey);}

	/**
	 * Returns the element index and master
	 * coodinates in grid/element space corresponding
	 * to the provided cartesian coordinate
	 */
	bool inter_inv_coords(const trafo_base &ct,double x,double y,
			int &ei,int &ej,double &xi,double &eta) const {

		// get intermediate coords
		double x_i,y_i;
		ct.x2xi(x,y,x_i,y_i);

		// if x or y is out of domain, return false
		if (x_i<ax || x_i>cx || y_i<ay || y_i>by) return false;

		// relate local coords to dist from "left" edge
		xi =  x_i > bx ? mi + (x_i-bx)/ex2 : (x_i-ax)/ex;
		eta = (y_i-ay)/ey;

		// integer part is index: store and subtract
		ei = xi>0  ? (int) xi  : (((int) xi)  - 1);
		xi-=ei;

		ej = eta>0 ? (int) eta : (((int) eta) - 1);
		eta-=ej;

		// the point was in the domain
		return true;
	}

	/**
	 * Returns element index and local / master coordinates
	 * at a given cartesian coordinate
	 */
	bool inv_coords(double x,double y,
			int &ei,int &ej,double &xi,double &eta) const {

		// if x or y is out of domain, return false
		if (x<ax || x>cx || y<ay || y>by) return false;

		// relate local coords to dist from "left" edge
		xi =  x > bx ? mi + (x-bx)/ex2 : (x-ax)/ex;
		eta = (y-ay)/ey;

		// integer part is index: store and subtract
		ei = xi>0  ? (int) xi  : (((int) xi)  - 1);
		xi-=ei;

		ej = eta>0 ? (int) eta : (((int) eta) - 1);
		eta-=ej;

		// the point was in the domain
		return true;
	}
};

/**
 * structure representing two regular cartesian grids glued
 * together at an intermediate x-value
 */
struct xsplit_grid_new: public cartesian_grid {

	/** index of the grid split */
	int mi1,mi2,mi3;
	/** "right"-most horizontal cartesian coord */
	double cx,dx,fx;
	/** element width in the second half of the split grid */
	double ex2,ex3,ex4;

	/**
	 * constructor
	 *
	 * @param[in] m1 number of columns of elements in the first (x) grid
	 * @param[in] m2 number of columns of elements in the second (x) grid
	 * @param[in] n  number of rows of elements in the (y) grid
	 * @param[in] ax left edge coordinate
	 * @param[in] bx intermediate coordinate
	 * @param[in] cx right edge coordinate
	 * @param[in] ay bottom edge coordinate
	 * @param[in] by top edge coordinate
	 */
	xsplit_grid_new(int m1,int m2,int m3,int m4,int n,
		double ax,double bx,double cx_,double dx_,double fx_,
		double ay,double by):
		cartesian_grid(m1,n,ax,bx,ay,by),
		mi1(m1),mi2(m2+m1),mi3(m1+m2+m3),
		cx(cx_),dx(dx_),fx(fx_),
		ex2((cx-bx)/m2),ex3((dx-cx)/m3),ex4((fx-dx)/m4) {}

	/** helper function returning the left x-coord of element i */
	double calc_ax(int i) {
		return ax + (
			i<mi1 ? (i*ex) : (
			i<mi2 ? (mi1*ex +   (i-mi1)*ex2) : (
			i<mi3 ? (mi1*ex + (mi2-mi1)*ex2 +   (i-mi2)*ex3) : (
				    (mi1*ex + (mi2-mi1)*ex2 + (mi3-mi2)*ex3 + (i-mi3)*ex4)
		))));
	}
	/** helper function returning the width of element i */
	double calc_ex(int i) {
		return i<mi1 ? ex : (
			   i<mi2 ? ex2 : (
			   i<mi3 ? ex3 : ex4
		));
	}

	/** returns the linear geom for element e */
	linear_geom element_geometry(element &e) {
		return linear_geom(e,calc_ax(e.i),ay+e.j*ey,calc_ex(e.i),ey);}
	/** allocates the linear geom for element e */
	trafo_base* alloc_element_geometry(element &e) {
		return new linear_geom(e,calc_ax(e.i),ay+e.j*ey,calc_ex(e.i),ey);}

	/**
	 * Returns the element index and master
	 * coodinates in grid/element space corresponding
	 * to the provided cartesian coordinate
	 */
	bool inter_inv_coords(const trafo_base &ct,double x,double y,
			int &ei,int &ej,double &xi,double &eta) const {

		// get intermediate coords
		double x_i,y_i;
		ct.x2xi(x,y,x_i,y_i);

		// if x or y is out of domain, return false
		if (x_i<ax || x_i>fx || y_i<ay || y_i>by) return false;

		// relate local coords to dist from "left" edge
		xi =  (x_i > dx) ? (mi3 + (x_i-dx)/ex4) : (
			  (x_i > cx) ? (mi2 + (x_i-cx)/ex3) : (
			  (x_i > bx) ? (mi1 + (x_i-bx)/ex2) : (
				                 (x_i-ax)/ex
		)));
		eta = (y_i-ay)/ey;

		// integer part is index: store and subtract
		ei = xi>0  ? (int) xi  : (((int) xi)  - 1);
		xi-=ei;

		ej = eta>0 ? (int) eta : (((int) eta) - 1);
		eta-=ej;

		// the point was in the domain
		return true;
	}

	/**
	 * Returns element index and local / master coordinates
	 * at a given cartesian coordinate
	 */
	bool inv_coords(double x,double y,
			int &ei,int &ej,double &xi,double &eta) const {

		// if x or y is out of domain, return false
		if (x<ax || x>fx || y<ay || y>by) return false;

		// relate local coords to dist from "left" edge
		xi =  x > dx ? (mi3 + (x-dx)/ex4) : (
			  x > cx ? (mi2 + (x-cx)/ex3) : (
		      x > bx ? (mi1 + (x-bx)/ex2) : (
				              (x-ax)/ex
		)));
		eta = (y-ay)/ey;

		// integer part is index: store and subtract
		ei = xi>0  ? (int) xi  : (((int) xi)  - 1);
		xi-=ei;

		ej = eta>0 ? (int) eta : (((int) eta) - 1);
		eta-=ej;

		// the point was in the domain
		return true;
	}
};

/**
 * struct representing a half-plane in polar coordinates,
 * inheriting from xsplit_grid xi x eta, with xi in (-xi1,0,xi2)
 * and eta in (0,pi), where xi1 and xi2 cooresponde to the surface
 * of two circles
 */
struct polar_half_plane: public cartesian_grid, public trafo_base {

	/** bipolar length scale */
	double a,ai;
	double R_max;
	double xi_max;

	/**
	 * Constructs a bipolar half-plane with (n x m1) points in the second
	 * quadrant and (n x m2) points in the first quadrant, such that the
	 * left- and right-most values correspond to the surfaces of two spheres
	 * with the provided radii
	 *
	 * @param[in] m1 number of element columns left of the y=0 line
	 * @param[in] m2 number of element columns right of the y=0 line
	 * @param[in] n number of element rows
	 * @param[in] d distance between the two spheres
	 * @param[in] R1 radius of the left sphere
	 * @param[in] R2 radius of the right sphere
	 */
	polar_half_plane(int m,int n,double a_,double R_max_):
		cartesian_grid(m,n,0,log(R_max_/a_),-M_PI,M_PI), // orig 0->PI
		a(a_),ai(1./a),R_max(R_max_),xi_max(log(R_max/a)) {}

	/** get real coordinates (x,y) at master coordinates (xi,eta) */
	void xi2x(double xi,double eta,double &x,double &y) const {
		double c = cos(eta);
		double s = sin(eta);

		x = a*exp(xi)*c; y = a*exp(xi)*s;}
	/** get determinant of the Jacobean at master coords (xi,eta) */
	double Jdet(double xi,double eta) const {return a*a*exp(2*xi);}
	/** get inverse Jacobean at (xi,eta) */
	void Jinv(double xi,double eta,
		 double &xi_x,double &et_x,double &xi_y,double &et_y) const {

		double c = cos(eta);
		double s = sin(eta);

		double fac = ai * exp(-xi);

		xi_x =  c * fac;
		et_x = -s * fac;
		xi_y =  s * fac;
		et_y =  c * fac;
	}
	void Jinv2(double xi,double et,
		double &xi_xx, double &xi_xy, double &xi_yy,
		double &et_xx, double &et_xy, double &et_yy) const {
		fprintf(stderr,"Jinv2 not implemented for circle\n");
		exit(1);
	}
	/** get bipolar coords at a cartesian coord */
	void x2xi(double x,double y,double &xi,double &eta) const {

		// radial dist
		double r = sqrt(x*x+y*y);
		
		xi  = log(r * ai);
		eta = atan2(y,x);
		if(eta < -M_PI) eta+=2*M_PI;
	}

	/**
	 * overload the geometry allocation using the square [xi,eta]
	 * grid as an intermediate to the bipolars
	 */
	trafo_base* alloc_element_geometry(element &e) {return coord_trafo(e,*this);}

	/**
	 * Returns the element index and master
	 * coodinates in grid/element space corresponding
	 * to the provided cartesian coordinate
	 */
	bool inv_coords(double x,double y,
			int &ei,int &ej,double &xi,double &eta) const {
		return inter_inv_coords(*this,x,y,ei,ej,xi,eta);
	}
};

/**
 * struct representing a half-plane in bipolar coordinates,
 * inheriting from xsplit_grid xi x eta, with xi in (-xi1,0,xi2)
 * and eta in (0,pi), where xi1 and xi2 cooresponde to the surface
 * of two circles
 */
struct bipolar_half_plane: public xsplit_grid_new, public trafo_base {

	/** bipolar length scale */
	double a;
	/** xi-coord of left,right sphere */
	double xil,xir;
	/** x-coord of left,right sphere */
	double xl,xr;
	/** helper function to get scale factor */
	static double drr2a(double d,double r1,double r2);
	/** helper function to get left-sphere xi */
	static double drr2xi1(double d,double r1,double r2);
	/** helper function to get right-sphere xi */
	static double drr2xi2(double d,double r1,double r2);
	/** returns the scale factor h = |dx|/|dxi| at xi,et */
	double h(double xi,double eta) const {return 1/(cosh(xi)-cos(eta));}

	/**
	 * Constructs a bipolar half-plane with (n x m1) points in the second
	 * quadrant and (n x m2) points in the first quadrant, such that the
	 * left- and right-most values correspond to the surfaces of two spheres
	 * with the provided radii
	 *
	 * @param[in] m1 number of element columns left of the y=0 line
	 * @param[in] m2 number of element columns right of the y=0 line
	 * @param[in] n number of element rows
	 * @param[in] d distance between the two spheres
	 * @param[in] R1 radius of the left sphere
	 * @param[in] R2 radius of the right sphere
	 */
	bipolar_half_plane(int m1,int m2,int m3,int m4,int n,double d,double R1,double R2,double bl_wid):

		xsplit_grid_new(m1,m2,m3,m4,n,

		drr2xi1(d,R1,R2),
		(
			(10 * bl_wid * drr2a(d,R1,R2) / R1 / R1) < 0.5*fabs(drr2xi1(d,R1,R2)) ?
			drr2xi1(d,R1,R2) + 10*bl_wid*drr2a(d,R1,R2)/R1/R1                     :
			0.5 * drr2xi1(d,R1,R2)
		),
		0,
		(
			(10 * bl_wid * drr2a(d,R1,R2) / R2 / R2) < 0.5*fabs(drr2xi2(d,R1,R2)) ?
			drr2xi2(d,R1,R2) - 10*bl_wid*drr2a(d,R1,R2)/R2/R2                     :
			0.5 * drr2xi2(d,R1,R2)
		),
		drr2xi2(d,R1,R2),
		-M_PI,M_PI),a(drr2a(d,R1,R2)),
		xil(drr2xi1(d,R1,R2)),xir(drr2xi2(d,R1,R2)),
		xl(a/tanh(xil)),xr(a/tanh(xir)) {}

	/** get real coordinates (x,y) at master coordinates (xi,eta) */
	void xi2x(double xi,double eta,double &x,double &y) const {
		x = a*h(xi,eta)*sinh(xi);y = a*h(xi,eta)*sin(eta);}
	/** get determinant of the Jacobean at master coords (xi,eta) */
	double Jdet(double xi,double eta) const
		{double hh=h(xi,eta); return a*a*hh*hh;}
	/** get inverse Jacobean at (xi,eta) */
	void Jinv(double xi,double eta,
		 double &xi_x,double &et_x,double &xi_y,double &et_y) const {
		et_y = -(xi_x = (1-cos(eta)*cosh(xi))/a);
		xi_y =   et_x =   -sin(eta)*sinh(xi) /a;

	}
	void Jinv2(double xi,double et,
		double &xi_xx, double &xi_xy, double &xi_yy,
		double &et_xx, double &et_xy, double &et_yy) const {

		double ia_sq = 1./(a*a);
		
		xi_xx = ( cosh(xi)*cos(2*et) - cos(et) ) * sinh(xi) * ia_sq;
		et_xx = ( 0.5*cosh(2*xi)*sin(2*et) - cosh(xi)*sin(et) ) * ia_sq;

		xi_xy = ( 0.5*cosh(2*xi)*sin(2*et) - cosh(xi)*sin(et) ) * ia_sq;
		et_xy = (-cosh(xi)*cos(2*et) + cos(et) ) * sinh(xi) * ia_sq;

		xi_yy = (-cosh(xi)*cos(2*et) + cos(et) ) * sinh(xi) * ia_sq;
		et_yy = (-0.5*cosh(2*xi)*sin(2*et) + cosh(xi)*sin(et) ) * ia_sq;
	}
	void Jinv_det(double xi,double eta,
		double &xi_x_xi,double &et_x_xi,double &xi_y_xi,double &et_y_xi,
		double &xi_x_et,double &et_x_et,double &xi_y_et,double &et_y_et) const{

		xi_x_xi = -cos(eta)*sinh(xi)/a;
		et_x_xi = -sin(eta)*cosh(xi)/a;
		xi_y_xi = -sin(eta)*cosh(xi)/a;
		et_y_xi =  cos(eta)*sinh(xi)/a;

		xi_x_et =  sin(eta)*cosh(xi)/a;
		et_x_et = -cos(eta)*sinh(xi)/a;
		xi_y_et = -cos(eta)*sinh(xi)/a;
		et_y_et = -sin(eta)*cosh(xi)/a;
	}
	/** get bipolar coords at a cartesian coord */
	void x2xi(double x,double y,double &xi,double &eta) const {

		// squared distances -> xi
		double dsq_l = (x+a)*(x+a) + y*y;
		double dsq_r = (x-a)*(x-a) + y*y;
		xi = 0.5*log(dsq_l/dsq_r);

		// dist ratios -> eta
		double v1 = fabs(2*a*y); // add sign back in
		double v2 = a*a - x*x - y*y;

		// numerator/denominator for arctan
		double n = v1, d = v2 + sqrt(v2*v2+v1*v1);
		eta = M_PI - 2*atan2(n,d);
		if (y<0) eta *= -1;
	}

	/**
	 * overload the geometry allocation using the square [xi,eta]
	 * grid as an intermediate to the bipolars
	 */
	trafo_base* alloc_element_geometry(element &e) {return coord_trafo(e,*this);}

	/**
	 * Returns the element index and master
	 * coodinates in grid/element space corresponding
	 * to the provided cartesian coordinate
	 */
	bool inv_coords(double x,double y,
			int &ei,int &ej,double &xi,double &eta) const {
		return inter_inv_coords(*this,x,y,ei,ej,xi,eta);
	}
};

/**
 * Represents a 3D axisymmetric domain corresponding to a 2D
 * half-plane rotated about one of the cartesian axes
 */
struct cylindrical: public domain {

	/** either rotated about x- or y-axis */
	enum axis {X,Y};
	/** axis of rotation */
	axis ax;
	/** half-plane domain */
	domain *dm;

	/**
	 * constructor takes in an axis of rotation and
	 * the 2D half-plane domain with which to rotate
	 */
	cylindrical(axis a_,domain *dm_):ax(a_),dm(dm_){}
	/** default is to delete rotated domain */
	~cylindrical() {if (dm!=NULL) delete dm;}

	/** trafo representing the revolved domain */
	struct revolve: public trafo_base {

		/** axis of rotation */
		axis ax;
		/** rotated trafo */
		trafo_base *t;

		/**
		 * constructor takes in an axis of rotation
		 * and a 2D trafo rule
		 */
		revolve(axis a_,trafo_base *t_):ax(a_),t(t_) {}
		/** delete trafo when done */
		~revolve() {if (t!=NULL) delete t;}

		/** get real coordinates (x,y) at master coordinates (xi,eta) */
		void xi2x(double xi,double eta,double &x,double &y) const {
			t->xi2x(xi,eta,x,y);}
		/**
		 * get determinant of the Jacobean at master coords (xi,eta)
		 *
		 * this is the only change due to the rotation: the determinant
		 * is multiplied by 2 Pi D, where D is the distance to the axis
		 * of rotation
		 */
		double Jdet(double xi,double eta) const {
			double x,y;
			t->xi2x(xi,eta,x,y);
			// PI -> 2PI
			double v=M_PI*t->Jdet(xi,eta) * fabs(ax==X?y:x);
			return v;
		}
		/** get inverse Jacobean at (xi,eta) */
		void Jinv(double xi,double eta,
			 double &xi_x,double &et_x,double &xi_y,double &et_y) const {
			t->Jinv(xi,eta,xi_x,et_x,xi_y,et_y);
		}
		void Jinv2(double xi,double et,
			double &xi_xx, double &xi_xy, double &xi_yy,
			double &et_xx, double &et_xy, double &et_yy) const {
			t->Jinv2(xi,et,xi_xx,xi_xy,xi_yy,et_xx,et_xy,et_yy);
		}
		void Jinv_det(double xi,double eta,
			double &c1,double &c2,double &c3,double &c4,
			double &c5,double &c6,double &c7,double &c8) const {
			t->Jinv_det(xi,eta,c1,c2,c3,c4,c5,c6,c7,c8);
		}

	};

	/**
	 * to allocate element geometry, first allocate underlying geom. it
	 * will be deleted when the revolve struct's destructor is called
	 */
	trafo_base* alloc_element_geometry(element &e) {
		return new revolve(ax,dm->alloc_element_geometry(e));
	}

	/**
	 * use revolved domain's global -> elem-wise coord trafo
	 */
	bool inv_coords(double x,double y,int &ei,int &ej,
			double &xi,double &eta) const {return dm->inv_coords(x,y,ei,ej,xi,eta);}
};

/**
 * Represents a 3D bispherical domain, corresponding to a
 * bipolar half-plane rotated about the x-axis
 */
struct bispherical: public cylindrical {

	/**
	 * Constructs a bispherical domain with the provided grid layout
	 * of quad elements in bipsherical space, such that the edges in
	 * xi-space correspond to the surfaces of two spheres of radii
	 * R1 and R2 whose centers are separated by a distance d
	 */
	bispherical(int m1,int m2,int m3,int m4,int n,double d,double R1,double R2,double i_bl=0):
		cylindrical(cylindrical::X,new bipolar_half_plane(m1,m2,m3,m4,n,d,R1,R2,
					i_bl > 0 ? (1./i_bl) : std::numeric_limits<double>::infinity()
					)) {}
};

/**
 * Represents a 3D bispherical domain, corresponding to a
 * bipolar half-plane rotated about the x-axis
 */
struct spherical: public cylindrical {

	/**
	 * Constructs a bispherical domain with the provided grid layout
	 * of quad elements in bipsherical space, such that the edges in
	 * xi-space correspond to the surfaces of two spheres of radii
	 * R1 and R2 whose centers are separated by a distance d
	 */
	spherical(int m,int n,double a,double R):
		cylindrical(cylindrical::X,new polar_half_plane(m,n,a,R)) {}
};

/**
 * Represents a 3D cylinder domain, corresponding to a
 * rectangular half-plane rotated about the x-axis
 */
struct cylinder: public cylindrical {

	/**
	 * Constructs a cylindrical domain with the provided grid layout
	 * of quad elements in cylindrical space consisting of a rotated
	 * rectangle of the provided dimensions
	 */
	cylinder(int m,int n,double Lx,double r):
		cylindrical(cylindrical::X,new cartesian_grid(m,n,-Lx/2,Lx/2,-r/2,r/2)) {}
};

/**
 * Represents a 2D wedge in polar space
 */
struct wedge : public cartesian_grid,public trafo_base {

	/**
	 * Constructs a wedge in polar space from the
	 * provided range of r and theta
	 */
	wedge(int m,int n,double r,double R,
			double th,double Th):cartesian_grid(m,n,r,R,th,Th) {}

	/** get real coordinates (x,y) at master coordinates (xi,eta) */
	void xi2x(double rr,double tt,double &x,double &y) const {
		x = rr * cos(tt);
		y = rr * sin(tt);
	}
	/** get determinant of the Jacobean at master coords (xi,eta) */
	double Jdet(double rr,double tt) const {return rr;}
	/** get inverse Jacobean at (xi,eta) */
	void Jinv(double rr,double tt,
		 double &rr_x,double &tt_x,double &rr_y,double &tt_y) const {
		
		rr_x =  cos(tt);
		rr_y =  sin(tt);
		tt_x = -sin(tt)/fabs(rr);
		tt_y =  cos(tt)/fabs(rr);
	}
	void Jinv2(double xi,double et,
		double &xi_xx, double &xi_xy, double &xi_yy,
		double &et_xx, double &et_xy, double &et_yy) const {
		fprintf(stderr,"Jinv2 not implemented for wedge\n");
		exit(1);
	}

	/** get element width cartesian_grid::coord_trafo */
	trafo_base* alloc_element_geometry(element &e) {return coord_trafo(e,*this);}
};

#endif
