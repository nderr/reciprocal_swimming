#include "fem.hh"

// elements

/**
 * helper function which sets nodal global indices
 * to reflect the values of the local indices
 */
void element::iterator::set_global() {
	// global node number can be formed out of
	// element coords and local node coords
	n.gi = ai + n.li;
	n.gj = aj + n.lj;
	check_index();

	n.g  = n.gi + um()*n.gj;
}

/**
 * increments the node in row-major fashion, by
 * incrementing the column. if the column index exits
 * the row, it is set to zero and the row is incremented
 */
element::iterator& element::iterator::operator++() {

	// local and column
	n.l++;n.li++;

	// check for row
	if (n.li > p) {n.li=0; n.lj++;}

	// adjust global index and return
	set_global();
	return *this;
}


/**
 * Fills a provided vector with all neighboring elements in the grid
 *
 * @param[out] es a vector with all neighboring elements in the grid
 */
void element::get_neighbors(std::vector<element> &es) {

	// loop through [-1,1] x [-1,1] grid
	for (int je=-1;je<2;je++) for (int ie=-1;ie<2;ie++)
			if (!(ie==0&&je==0)) {

		// add to vector if in the grid
		element n = neigh(ie,je);
		if (n.in_grid()) es.push_back(n);
	}
}

/**
 * Returns a vector full of the shared nodes between this
 * and a provided element. The vector constists of node pairs,
 * such that the first in the pair this element's node and
 * the second is the other's.
 *
 * @param[in] other the element being searched
 * @param[out] shared a vector of all nodes shared with other
 */
void element::get_shared_nodes(element other,std::vector<node_pair> &shared) {
	// loop through each element's nodes
	for (node me : *this) for (node you : other)

		// store if nodes have the same global index
		if (me.g == you.g) shared.emplace_back(me,you);
}

// grid


/**
 * increments the column index. if the column index
 * exits the row, sets column index to zero and increments
 * the row index
 */
grid::iterator& grid::iterator::operator++() {
	e.i++;
	if (e.i >= g->m) {e.i=0; e.j++;}
	return *this;
}

/**
 * decrements the column index. if it becomes negative,
 * sets column index to the last in the row and decrements
 * the row index
 */
grid::iterator& grid::iterator::operator--() {
	e.i--;
	if (e.i < 0) {e.i=g->m-1; e.j--;}
	return *this;
}

/**
 * addition overload, allowing pointer-type behavior,
 * e.g. *(it+3) returns the element 3 after the one
 * pointed to by an iterator it
 */
grid::iterator grid::iterator::operator+(int n) const {

	// extract sign and absolute value of addend
	int sgn = n>=0?1:-1;
	n = abs(n);

	// use increments if possible since they're cheaper
	if (n==1) {
		iterator other = *this;
		return sgn>0 ? (++other) : (--other);
	}

	// calculate new col/row indices
	int new_i = e.i + sgn * (n%g->m);
	int new_j = e.j + sgn * (n/g->m);

	// adjust until column index is in range, then return
	while (new_i >= g->m) {new_i -= g->m; new_j++;}
	while (new_i < 0)     {new_i += g->m; new_j--;}
	return iterator(*g,new_i,new_j);
}

/** add/equals operator procees much the same as addition */
grid::iterator grid::iterator::operator+=(int n) {

	// extract sign/absolute value
	int sgn = n>=0?1:-1;
	n = abs(n);

	// use increment / decrement if possible since cheaper
	if (n==1) return sgn==1?(++(*this)):(--(*this));

	// calculate new col/row
	e.i += sgn * (n%g->m);
	e.j += sgn * (n/g->m);

	// adjust column until in range
	while (e.i >= g->m) {e.i -= g->m; e.j++;}
	while (e.i  < 0)    {e.i += g->m; e.j--;}
	return *this;
}

/** non-member addition for grid iterators */
grid::iterator operator+(int n,const grid::iterator it) {return it+n;}
/** non-member subtraction for grid iterators */
grid::iterator operator-(int n,const grid::iterator it) {return it-n;}

/**
 * increments the column index. if the column index
 * exits the row, sets column index to zero and increments
 * the row index
 */
gnode_grid::iterator& gnode_grid::iterator::operator++() {
	n.gi++;
	if (n.gi >= g->m) {n.gi=0; n.gj++;}
	return *this;
}

/**
 * decrements the column index. if it becomes negative,
 * sets column index to the last in the row and decrements
 * the row index
 */
gnode_grid::iterator& gnode_grid::iterator::operator--() {
	n.gi--;
	if (n.gi < 0) {n.gi=g->m-1; n.gj--;}
	return *this;
}

/**
 * addition overload, allowing pointer-type behavior,
 * e.g. *(it+3) returns the element 3 after the one
 * pointed to by an iterator it
 */
gnode_grid::iterator gnode_grid::iterator::operator+(int d) const {

	// extract sign and absolute value of addend
	int sgn = d>=0?1:-1;
	d = abs(d);

	// use increments if possible since they're cheaper
	if (d==1) {
		iterator other = *this;
		return sgn>0 ? (++other) : (--other);
	}

	// calculate new col/row indices
	int new_i = n.gi + sgn * (d%g->m);
	int new_j = n.gj + sgn * (d/g->m);

	// adjust until column index is in range, then return
	while (new_i >= g->m) {new_i -= g->m; new_j++;}
	while (new_i < 0)     {new_i += g->m; new_j--;}
	return iterator(*g,new_i,new_j);
}

/** add/equals operator procees much the same as addition */
gnode_grid::iterator gnode_grid::iterator::operator+=(int d) {

	// extract sign/absolute value
	int sgn = d>=0?1:-1;
	d = abs(d);

	// use increment / decrement if possible since cheaper
	if (d==1) return sgn==1?(++(*this)):(--(*this));

	// calculate new col/row
	n.gi += sgn * (d%g->m);
	n.gj += sgn * (d/g->m);

	// adjust column until in range
	while (n.gi >= g->m) {n.gi -= g->m; n.gj++;}
	while (n.gi  < 0)    {n.gi += g->m; n.gj--;}
	return *this;
}

/** non-member addition for grid iterators */
gnode_grid::iterator operator+(int d,const gnode_grid::iterator it) {return it+d;}
/** non-member subtraction for grid iterators */
gnode_grid::iterator operator-(int d,const gnode_grid::iterator it) {return it-d;}

// cartesian_grid

/** get real coordinates (x,y) at master coordinates (xi,eta) */
void cartesian_grid::trafo_geom::xi2x(double xi,double et,double &x,double &y) const {

	// grab intermediate coords, then plug them in
	double x_i,y_i;
	p.xi2x(xi,et,x_i,y_i);
	ct.xi2x(x_i,y_i,x,y);
}

/** get determinant of the Jacobean at master coords (xi,eta) */
double cartesian_grid::trafo_geom::Jdet(double xi,double et) const {

	// grab intermediate coords, then plug them in
	double x_i,y_i;
	p.xi2x(xi,et,x_i,y_i);
	return p.Jdet(xi,et) * ct.Jdet(x_i,y_i);
}

/** get inverse Jacobean at (xi,eta) */
void cartesian_grid::trafo_geom::Jinv(double xi,double et,double &xi_x,
		double &et_x,double &xi_y,double &et_y) const {

	// only need two of the intermediate Jinv terms
	double xi_x_i = 1/p.dx, et_y_i = 1/p.dy;

	// grab intermediate coords, then plug them in
	double x_i,y_i;
	p.xi2x(xi,et,x_i,y_i);

	// grab inverse jacobean...
	double x_i_x,y_i_x,x_i_y,y_i_y;
	ct.Jinv(x_i,y_i,x_i_x,y_i_x,x_i_y,y_i_y);

	// ...chain rule
	xi_x = xi_x_i*x_i_x;
	et_x = et_y_i*y_i_x;
	xi_y = xi_x_i*x_i_y;
	et_y = et_y_i*y_i_y;

//	printf("xi_x_i: %g, x_i_x: %g, xi_x: %g\n",xi_x_i,x_i_x,xi_x);
//	printf("et_y_i: %g, y_i_x: %g, et_x: %g\n",et_y_i,y_i_x,et_x);
}

/** get inverse Jacobean at (xi,eta) */
void cartesian_grid::trafo_geom::Jinv2(double xi,double et,
	double &xi_xx, double &xi_xy, double &xi_yy,
	double &et_xx, double &et_xy, double &et_yy) const {

	// only need two of the intermediate Jinv terms
	double xi_x_i = 1/p.dx, et_y_i = 1/p.dy;

	// grab intermediate coords, then plug them in
	double x_i,y_i;
	p.xi2x(xi,et,x_i,y_i);

	// grab inverse jacobean...
	double x_i_xx, y_i_xx, x_i_xy, y_i_xy, x_i_yy, y_i_yy;
	ct.Jinv2(x_i,y_i,x_i_xx,x_i_xy,x_i_yy,y_i_xx,y_i_xy,y_i_yy);

	// ...chain rule
	xi_xx = xi_x_i * x_i_xx; et_xx = et_y_i * y_i_xx;
	xi_xy = xi_x_i * x_i_xy; et_xy = et_y_i * y_i_xy;
	xi_yy = xi_x_i * x_i_yy; et_yy = et_y_i * y_i_yy;

//	printf("xi_x_i: %g, x_i_x: %g, xi_x: %g\n",xi_x_i,x_i_x,xi_x);
//	printf("et_y_i: %g, y_i_x: %g, et_x: %g\n",et_y_i,y_i_x,et_x);
}

/** get inverse Jacobean at (xi,eta) */
void cartesian_grid::trafo_geom::Jinv_det(double xi,double et,
		double &xi_x_xi,double &et_x_xi,double &xi_y_xi,double &et_y_xi,
		double &xi_x_et,double &et_x_et,double &xi_y_et,double &et_y_et) const{

	// only need two of the intermediate Jinv terms
	double xi_x_i = 1/p.dx, et_y_i = 1/p.dy;

	// grab intermediate coords, then plug them in
	double x_i,y_i;
	p.xi2x(xi,et,x_i,y_i);

	// grab inverse jacobean dets...
	//double x_i_x,y_i_x,x_i_y,y_i_y;
	double x_i_x_xi,y_i_x_xi,x_i_y_xi,y_i_y_xi,x_i_x_et,y_i_x_et,x_i_y_et,y_i_y_et;
	ct.Jinv_det(x_i,y_i,x_i_x_xi,y_i_x_xi,x_i_y_xi,y_i_y_xi,
			x_i_x_et,y_i_x_et,x_i_y_et,y_i_y_et);

	// ...chain rule
	//xi_x = xi_x_i*x_i_x;
	//et_x = et_y_i*y_i_x;
	//xi_y = xi_x_i*x_i_y;
	//et_y = et_y_i*y_i_y;
	xi_x_xi = xi_x_i * x_i_x_xi;
	et_x_xi = et_y_i * y_i_x_xi;
	xi_y_xi = xi_x_i * x_i_y_xi;
	et_y_xi = et_y_i * y_i_y_xi;

	xi_x_et = xi_x_i * x_i_x_et;
	et_x_et = et_y_i * y_i_x_et;
	xi_y_et = xi_x_i * x_i_y_et;
	et_y_et = et_y_i * y_i_y_et;
}

/**
 * Returns the length scale corresponding to the bipolar
 * coordinate system such that two spheres of radii r1
 * and r2 have a distance d between their centers
 *
 * @param[in] d distance between spheres
 * @param[in] r1 radius of left sphere
 * @param[in] r2 radius of right sphere
 * @return the scale factor of the coordinate system
 */
double bipolar_half_plane::drr2a(double d,double r1,double r2) {

	// calculate intermediate values
	double v1=(r1+r2)/d, v2=2*r1*r2/(d*d), v=1-v1*v1+v2;

	// plug back in
	return 0.5*d*sqrt(v*v-v2*v2);
}

/**
 * Returns the xi-value corresponding to the surface
 * of the left sphere
 *
 * @param[in] d distance between spheres
 * @param[in] r1 radius of left sphere
 * @param[in] r2 radius of right sphere
 * @return the scale factor of the coordinate system
 */
double bipolar_half_plane::drr2xi1(double d,double r1,double r2) {

	// calculate scale factor
	double a = drr2a(d,r1,r2);

	// xi is a function of scale factor and radius
	return -asinh(a/r1);
}

/**
 * Returns the xi-value corresponding to the surface
 * of the right sphere
 *
 * @param[in] d distance between spheres
 * @param[in] r1 radius of left sphere
 * @param[in] r2 radius of right sphere
 * @return the scale factor of the coordinate system
 */
double bipolar_half_plane::drr2xi2(double d,double r1,double r2) {

	// calculate scale factor
	double a = drr2a(d,r1,r2);

	// xi is a function of scale factor and radius
	return asinh(a/r2);
}
