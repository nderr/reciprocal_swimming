#include "master.hh"

/**
 * Cleans up helper objects
 */
master::~master() {
	cleanup_shape_functions();
	cleanup_quad();
}

/**
 * Returns the local index of the node's corners
 *
 * e.g. p=2:
 *
 *     5
 *    /|
 *   3 4
 *  /  |
 * 0-1-2
 *
 * @param[out] crn pointer to set of integers to be set
 *   to the corners' node indices
 */
void master_tri::corners(int *crn) {
	(*crn++)=0;
	(*crn++)=p;
	(*crn++)=n-1;
}

/**
 * Returns the local index of the node's corners
 *
 * e.g. p=2:
 *
 * 6-7-8
 * | | |
 * 3-4-5
 * | | |
 * 0-1-2
 *
 * @param[out] crn pointer to set of integers to be set
 *   to the corners' node indices
 */
void master_quad::corners(int *crn) {
	(*crn++)=0;
	(*crn++)=p;
	(*crn++)=n-1;
	(*crn++)=n-1-p;
}

/**
 * Returns the local indices of one of the node's sides
 *
 * e.g. p=2:
 *
 *     5
 *    /|
 *   3 4
 *  /  |
 * 0-1-2
 *
 *  si | sid
 * ----|-----------
 *  0  | 0,1,2
 *  1  | 2,4,5
 *  2  | 0,3,5
 *
 * @param[in] si the index of the side in question
 * @param[out] crn pointer to set of integers to be set
 *   to the specified sides' node indices
 */
void master_tri::side(int si,int *sid) {
	int base;
	switch (si) {
		case 0:
			for(int i=0;i<=p;i++) (*sid++)=i;
			break;
		case 1:
			base=0;
			for(int i=0;i<=p;i++) {
				(*sid++)=base+(p-i);
				base+=p-i+1;
			}
			break;
		case 2:
			base=0;
			for(int i=0;i<=p;i++) {
				(*sid++)=base;
				base += p-i+1;
			}
			break;
	}
}

/**
 * Returns the local indices of one of the node's sides
 * 
 * e.g. p=2:
 * 
 * 6-7-8
 * | | |
 * 3-4-5
 * | | |
 * 0-1-2
 *
 *  si | sid
 * ----|-----------
 *  0  | 0,1,2
 *  1  | 2,5,8
 *  2  | 8,7,6
 *  3  | 6,3,0
 *
 * @param[in] si the index of the side in question
 * @param[out] crn pointer to set of integers to be set
 *   to the specified sides' node indices
 */
void master_quad::side(int si,int *sid) {
	switch (si) {
		case 0:
			for(int i=0;i<=p;i++) (*sid++)=i;
			break;
		case 1:
			for(int i=0;i<=p;i++) (*sid++)=i*(p+1)+p;
			break;
		case 2:
			for(int i=0;i<=p;i++) (*sid++)=n-1-i;
			break;
		case 3:
			for(int i=0;i<=p;i++) (*sid++)=(p-i)*(p+1);
			break;
	}
}

/**
 * Returns the coordinates of the ith local node in xi-eta space
 *
 * @param[in] i index of the requested node
 * @param[out] xi the xi-coordinate of the requested node
 * @param[out] eta the eta-coordiante of the requested node
 */
void master_tri::point(int i,double &xi,double &eta) {

	// nodes are linearly distributed throughout the element
	//
	// XXX one possible improvement: replace w Chebyshev
	for (int k=0,ii=0;k<=p;k++) {
		for (int j=0;j<=(p-k);j++,ii++) {
			if(ii==i) {
				xi=static_cast<double>(j)/p;
				eta=static_cast<double>(k)/p;
				return;
			}
		}
	}

	fputs("bad indexing!",stderr);
	exit(-1);
}

/**
 * Returns the coordinates of the ith local node in xi-eta space
 *
 * @param[in] i index of the requested node
 * @param[out] xi the xi-coordinate of the requested node
 * @param[out] eta the eta-coordiante of the requested node
 */
void master_quad::point(int i,double &xi,double &eta) {

	// nodes are linearly distributed throughout the element
	//
	// XXX one possible improvement: replace w Chebyshev
	int j,k;
	k = i/(p+1);
	j = i%(p+1);
	xi = static_cast<double>(j)/p;
	eta = static_cast<double>(k)/p;
}

/**
 * Calculates shape functions given polynomial set and
 * quadrature points. This function should be called after
 * the gauss_quad_2d object is allocated
 */
void master::setup_shape_functions() {

	// allocate new polynomial object
	// and set node number accordingly
	//poly_2d *polys = get_new_poly();
	polys = get_new_poly();
	n = polys->n;

	// allocate matrices

	// n x n (ith base polynomial,
	//        jth Lagrange polynomial)
	lagr = new double[n*(n+6*g->N)];

	// g x n (ith gauss quad point,
	//        jth basis polynomial)
	sh    = lagr+n*n;
	sh_xi = sh+n*g->N;
	sh_et = sh_xi+n*g->N;
	sh_L = sh_et+n*g->N;
	sh_L_xi = sh_L+n*g->N;
	sh_L_et = sh_L_xi+n*g->N;


	// get coefficients of Lagrange polynomials by setting up
	// matrix V_{ji} = p_i(x_j) [i.e. row: point, col: basis func]
	// and inverting in place so that coefficients C_{ij}
	// [i.e. row: function, col: coeff] satisfies 
	// C_{ij} p_j(x_k) = delta{ik} => C=V^{-T}.
	//
	// Note that we lagr corresponds to C^T, since it's the inverse,
	// not transpose inverse, of V.
	//
	// Since the shape functions are stored in g x n arrays M_{ij}=p_j(x_i),
	// where p_j is the jth basis polynomial, keeping c=C^T is convenient
	// since M_{ij} c_{jk} is a (g x n) array m_{ij} =L_j(x_i) where L_j is the
	// jth Lagrange polynomial
	//
	// TODO just store this latter matrix m_{ij}! Neither M nor c are used
	// anywhere else.
	double xi,et;
	for (int j=0;j<n;j++) {
		point(j,xi,et);
		for (int i=0;i<n;i++) {
			lagr[i+j*n] = polys->eval(i,2*xi-1,2*et-1);
		}
	}
	lapack::invert(n,lagr);

	// set up gauss quadrature weighting. the polynomial object
	// defines said functions on the [-1,1]^d space, so node values
	// in [0,1]^d space must be converted before befor evaluation.
	// Derivative values must be doubled.
	for (int i=0;i<g->N;i++) {
		g->point(i,xi,et);
		for (int j=0;j<n;j++) {
			polys->eval(j,2*xi-1,2*et-1,sh[j+i*n],sh_xi[j+i*n],sh_et[j+i*n]);
			sh_xi[j+i*n]*=2;
			sh_et[j+i*n]*=2;
		}
	}

	// set up shape functions w.r.t. Lagrange rather than basis polynomials
	for (int i=0;i<g->N;i++) {

		// at each gauss point i, look at a single node j
		for (int j=0;j<n;j++) {

			// grab references to Lagrange shape functions at this point
			double &v=sh_L[j+i*n];
			double &v_xi=sh_L_xi[j+i*n];
			double &v_et=sh_L_et[j+i*n];

			// zero them out
			v=v_xi=v_et=0;

			// now get contribution from each basis function
			for (int k=0;k<n;k++) {

				// get the contribution of the kth basis polynomial
				// to the jth Lagrange polynomial
				double L = lagr[j+k*n];

				// get the jth Lagrange poly at the ith gauss node
				// by summing up the contributions of each basis poly
				v    +=    sh[k+i*n]*L;
				v_xi += sh_xi[k+i*n]*L;
				v_et += sh_et[k+i*n]*L;
			}
		}
	}

	// clean up polynomial object
//	delete polys; XXX now deleted in cleanup_shape_functions
}

/**
 * Prints out mass matrix of unscaled master element
 * for verification purposes
 */
void master::verify() {

	// allocate matrix
	double *m = new double[n*n];
	for (int i=0;i<n*n;i++) m[i]=0;

	// grab each gauss node
	double x,y;
	for (int i=0;i<g->N;i++) {
		g->point(i,x,y);

		// construct mass matrix using quadrature for integrals
		for(int j=0;j<n;j++) for(int k=0;k<n;k++) {
			m[j+k*n]+=g->weight(i)*sh[j+i*n]*sh[k+i*n];
		}
	}

	// print matrix and clean up
	print_mat(m,n,n);
	delete[] m;
}

/**
 * Outputs the Stokes-style mass matrix for the mixed finite
 * element problem (v,p)
 *
 * TODO replace with stokes_master subclass
 */
/*
void master::stokes_mat(master *pmast,trafo *ct,double *m) {

	// check for polynomial order consistency
	if (ng!=pmast->ng) {
		printf("%d vs %d\n",ng,pmast->ng);
		fputs("inconsistent quadrature points\n",stderr);
		exit(-1);
	}

	// zero output
	for(int i=0;i<2*n*pmast->n;i++) m[i]=0;

	// preallocate variables
	double Ji_xi_x,Ji_xi_y,Ji_et_x,Ji_et_y;
	double sh_xi_j,sh_et_j;
	double dx_j,dy_j;
	double sh_k;
	double wght,det;
	double xi,et;
	for (int i=0;i<g->N;i++) {

		g->point(i,xi,et);
		wght = g->weight(i);
		det = ct->J(xi,et,Ji_xi_x,Ji_xi_y,Ji_et_x,Ji_et_y);

		for(int j=0;j<n;j++) {

			sh_xi_j = 0;
			sh_et_j = 0;
			for(int si=0;si<n;si++) {
				sh_xi_j += sh_xi[si+i*n]*lagr[j+si*n];
				sh_et_j += sh_et[si+i*n]*lagr[j+si*n];
			}

			dx_j = Ji_xi_x*sh_xi_j + Ji_et_x*sh_et_j;
			dy_j = Ji_xi_y*sh_xi_j + Ji_et_y*sh_et_j;

			for(int k=0;k<pmast->n;k++) {

				sh_k = 0;
				for(int si=0;si<pmast->n;si++) {
					sh_k += pmast->sh[si+i*pmast->n]
						* pmast->lagr[k+si*pmast->n];
				}

				m[k +   (2*j)*pmast->n] += wght*(dx_j*sh_k)*det;
				m[k + (2*j+1)*pmast->n] += wght*(dy_j*sh_k)*det;
			}
		}
	}
}
*/

/**
 * Output element boundaries to be plotted as lines in gnuplot
 *
 * @param[in] ct coordinate transform defining real element geometry
 * @param[in] out filestream to print to
 */
void master::gp_lines(const trafo_base &ct,double *dat,FILE *out) {

	// preallocate variables
	double xi,et,x,y;

	// pull in side indices
	int *sid = new int[p+1];

	// for each side, plot segments
	for(int i=0;i<=n_sides();i++) {
		side(i,sid);
		for(int j=0;j<=p;j++) {
			point(sid[j],xi,et);
			ct.xi2x(xi,et,x,y);

			fprintf(out,"%g %g %g\n",x,y,dat[sid[j]]);
		}
		fprintf(out,"\n\n");
	}

	// cleanup
	delete[] sid;
}

/**
 * Output scalar function for plotting within gnuplot
 *
 * TODO this is currently only compatible with quad elements. make general
 *
 * e.g. for p=1
 *
 * 2-3
 * | |
 * 0-1
 *
 * with data dat = (d0,d1,d2,d3)...
 *
 * output looks like
 *
 * SOF
 * x0 y0 d0
 * x1 y1 d1
 *
 * x2 y2 d2
 * x3 y3 d3
 *
 *
 * EOF
 *
 * @param[in] ct coordinate transform defining real element geometry
 * @param[in] dat coefficients of local nodes representing scalar functions
 *   to plot
 * @param[in] out filestream to print to
 */
void master::gp_fill(const trafo_base &ct,double *dat,FILE *out) {

	// pull in list of indices on the right side
	int *newlines = new int[p+1];
	side(1,newlines);

	// preallocate vars and loop over nodes
	double xi,eta,x,y;
	for(int i=0,nl=0;i<n;i++) {

		// grab coordinates and print out 
		point(i,xi,eta);
		ct.xi2x(xi,eta,x,y);
		fprintf(out,"%g %g %g\n",x,y,dat[i]);

		// drop a new line if we're at the end of a row of values
		if (i==newlines[nl]) {
			fprintf(out,"\n");
			nl++;
		}
	}

	// drop a final newline and cleanup
	fprintf(out,"\n");
	delete[] newlines;
}

void master::check_stokes_up(master &mu,master &mp) {
	if (mu.ng != mp.ng || mu.p != mp.p+1) {
		fprintf(stderr,"gauss points: u(%d), p(%d)\n",mu.ng,mp.ng);
		fprintf(stderr,"poly order  : u(%d), p(%d)\n",mu.p,mp.p);
		fputs("bad u vs match for Stokes problem\n",stderr);
		exit(-1);
	}
}

void master::stokes_bmats(master &mu,master &mp,const trafo_base &ct,double *bx,
		double *by,bool cyl,bool add) {

	check_stokes_up(mu,mp);
	
	if (!add) for(int i=0;i<mu.n*mp.n;i++) {bx[i]=by[i]=0;}

	// preallocate variables
	double wght,det;
	double xi,et;

	// u is rows
	double u_xi_j,u_et_j,u_j;
	double Ji_xi_x,Ji_et_x,Ji_xi_y,Ji_et_y;
	double u_x_j,u_y_j;

	// p is cols
	double p_k;
	double x,y;
	
	// loop through gauss points
	for (int i=0;i<mu.g->N;i++) {

		mu.g->point(i,xi,et);
		wght = mu.g->weight(i);
		det = ct.Jdet(xi,et);
		ct.Jinv(xi,et,Ji_xi_x,Ji_et_x,Ji_xi_y,Ji_et_y);
		ct.xi2x(xi,et,x,y);

		// loop through rows
		for (int j=0;j<mu.n;j++) {

			u_j    = mu.sh_L[j+i*mu.n];

			u_xi_j = mu.sh_L_xi[j+i*mu.n];
			u_et_j = mu.sh_L_et[j+i*mu.n];

			u_x_j = u_xi_j*Ji_xi_x + u_et_j*Ji_et_x;
			u_y_j = u_xi_j*Ji_xi_y + u_et_j*Ji_et_y;

			// loop through cols
			for (int k=0;k<mp.n;k++) {

				p_k = mp.sh_L[k+i*mp.n];

				// put in matrix
				bx[k + j*mp.n] += wght*u_x_j*p_k*det;
				by[k + j*mp.n] += wght*u_y_j*p_k*det;

				if (cyl) by[k + j*mp.n] += wght*u_j*p_k*det/y;
			}
		}
	}
}

/**
 * Calculates an element's mass matrix and stores it
 * in a provided array, provided a coordinate transform
 * 
 * @param[in] CT coordinate tranform object
 * @param[in] ct coordinate transform determining an
 *   element's real geometry
 * @param[out] m pointer to double array into which mass
 *   matrix is written
 */
void master::stiffness_matrix(const trafo_base &ct,double *m,bool add) {

	// zero output
	if (!add) for(int i=0;i<n*n;i++) m[i]=0;

	// preallocate variables
	double wght,det;
	double xi,et;
	double sh_j,sh_k;
	double sh_xi_j,sh_xi_k,sh_et_j,sh_et_k;
	double Ji_xi_x,Ji_xi_y,Ji_et_x,Ji_et_y;
	double dx_j,dx_k,dy_j,dy_k;

	// march though gauss nodes
	for (int i=0;i<g->N;i++) {

		// grab coordinate, weight and determinant
		g->point(i,xi,et);
		wght = g->weight(i);

		// get inverse Jacobean and determinatn
		det = ct.Jdet(xi,et);
		ct.Jinv(xi,et,Ji_xi_x,Ji_et_x,Ji_xi_y,Ji_et_y);

		double pre=wght*det;

//		printf("%g %g\n",wght,det);

		// march through element nodes
		for(int j=0;j<n;j++) {

			sh_j = sh_L[j+i*n];

			// get shape function derivative values for
			// the jth Lagrange polynomial at ith quad point
			sh_xi_j = sh_L_xi[j+i*n];
			sh_et_j = sh_L_et[j+i*n];

			// convert to real geometry derivatives
			// of jth Lagrange polynomial
			dx_j = Ji_xi_x*sh_xi_j + Ji_et_x*sh_et_j;
			dy_j = Ji_xi_y*sh_xi_j + Ji_et_y*sh_et_j;

			// march through element nodes
			for(int k=0;k<n;k++) {

				// get shape function derivative values for
				// the kth Lagrange polynomial
				sh_xi_k = sh_L_xi[k+i*n];
				sh_et_k = sh_L_et[k+i*n];
				
				// convert to real geometry derivatives of
				// kth Lagrange polynomial
				dx_k = Ji_xi_x*sh_xi_k + Ji_et_x*sh_et_k;
				dy_k = Ji_xi_y*sh_xi_k + Ji_et_y*sh_et_k;

				// record the contribution of j/k node combination
				// to the element mass matrix
				m[k+j*n] += pre*(dx_j*dx_k+dy_j*dy_k);
			}
		}
	}
}

/**
 * Calculates an element's mass matrix and stores it
 * in a provided array, provided a coordinate transform
 * 
 * @param[in] CT coordinate tranform object
 * @param[in] ct coordinate transform determining an
 *   element's real geometry
 * @param[out] m pointer to double array into which mass
 *   matrix is written
 */
void master::mass_matrix_cyl_y(const trafo_base &ct,double *m,bool add) {

	// zero output
	if (!add) for(int i=0;i<n*n;i++) m[i]=0;

	// preallocate variables
	double wght,det;
	double xi,et;
	double sh_j,sh_k;

	// march though gauss nodes
	for (int i=0;i<g->N;i++) {

		// grab coordinate, weight and determinant
		g->point(i,xi,et);
		wght = g->weight(i);

		// get inverse Jacobean and determinatn
		det = ct.Jdet(xi,et);

		double pre=wght*det;

		// weight mass matrix by 1/y^2
		double x,y;
		ct.xi2x(xi,et,x,y);
		pre /= (y*y);

		// march through element nodes
		for(int j=0;j<n;j++) {

			// get shape function values for the jth Lagrange polynomial at ith quad point
			sh_j = sh_L[j+i*n];

			// march through element nodes
			for(int k=0;k<n;k++) {

				// get shape function values for the kth Lagrange polynomial
				sh_k = sh_L[k+i*n];

				// record the contribution of j/k node combination
				// to the element mass matrix
				m[k+j*n] += pre*sh_j*sh_k;
			}
		}
	}
}

/**
 * Calculates an element's mass matrix and stores it
 * in a provided array, provided a coordinate transform
 * 
 * @param[in] CT coordinate tranform object
 * @param[in] ct coordinate transform determining an
 *   element's real geometry
 * @param[out] m pointer to double array into which mass
 *   matrix is written
 */
void master::mass_matrix(const trafo_base &ct,double *m,bool add) {

	// zero output
	if (!add) for(int i=0;i<n*n;i++) m[i]=0;

	// preallocate variables
	double wght,det;
	double xi,et;
	double sh_j,sh_k;

	// march though gauss nodes
	for (int i=0;i<g->N;i++) {

		// grab coordinate, weight and determinant
		g->point(i,xi,et);
		wght = g->weight(i);

		// get inverse Jacobean and determinatn
		det = ct.Jdet(xi,et);

		double pre=wght*det;

		// march through element nodes
		for(int j=0;j<n;j++) {

			// get shape function values for the jth Lagrange polynomial at ith quad point
			sh_j = sh_L[j+i*n];

			// march through element nodes
			for(int k=0;k<n;k++) {

				// get shape function values for the kth Lagrange polynomial
				sh_k = sh_L[k+i*n];

				// record the contribution of j/k node combination
				// to the element mass matrix
				m[k+j*n] += pre*sh_j*sh_k;
			}
		}
	}
}

/**
 * Calculates a source term S = int -f(x)v(x) dV = s_i v_i given a source
 * function f(x) and with v(x) = v_i phi_i(x), with phi_i the ith
 * source function.
 *
 * The integral above in x-space is here computed in xi-space, given
 * a coordinate transform object detailing the x <-> xi transformation
 *
 * TODO convert to template function taking functor for source
 *
 * @param[in] ct coordinate transform object defining real
 *   geometry of element
 * @param[in] sfun source function object
 * @param[out] s nx1 vector corresponding to integrated source term
 */
void master::source_term(const trafo_base &ct,const src_base &src,double *s) {

	// zero output vector
	for(int i=0;i<n;i++) s[i]=0;

	// preallocate variables
	double wght,det;
	double xi,et,x,y;
	double sh_pt,src_pt;

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

		double pre=wght*det*src_pt;

		// loop through vector elements
		for(int j=0;j<n;j++) {

			// get contribution of jth Lagrange function at
			// ith quadrature point
			sh_pt = sh_L[j+i*n];

			// add ith quadrature point's contribution
			// to jth vector element
			s[j] += pre*sh_pt;
		}
	}
}
