///////////////////////////////////////////////////////////////////////////////
// gp_data.cc                                                                //
//                                                                           //
// code for processing / outputing gnuplot graphics from a brinkman FEM      //
// solution as saved in a binary file                                        //
///////////////////////////////////////////////////////////////////////////////

#include "gp_data.hh"

/**
 * Returns the part of a complex number corresponding
 * to the provided C-string:
 *
 * r - real
 * i - imag
 * n - norm (magnitude in x-y plane)
 * a - argument
 *
 * @param[in] s string code of complex part
 * @return specified complex part
 */
gp_data::complex_part gp_data::read_part(const std::string &str) {

	if      (str=="r") return complex_part::real;
	else if (str=="i") return complex_part::imag;
	else if (str=="n") return complex_part::norm;
	else if (str=="a") return complex_part::arg;

	std::cerr << "bad part code " << str << std::endl;
	exit(-1);
}

/**
 * Returns the scalar field corresponding
 * to the provided C-string:
 *
 * ux - horizontal velocity
 * uy - vertical velocity
 * um - velocity magnitude
 * p  - pressure 
 * w  - vorticity
 * d  - divergence
 *
 * @param[in] s string code of scalar field
 * @return specified scalar
 */
gp_data::scalar gp_data::read_scalar(const std::string &str) {

	if      (str=="ux")    return scalar::vel_x;
	else if (str=="uy")    return scalar::vel_y;
	else if (str=="um")    return scalar::vel_mag;
	else if (str=="p")     return scalar::pres;
	else if (str=="w")     return scalar::vort;
	else if (str=="d")     return scalar::div;
	else if (str=="none")  return scalar::none;

	std::cerr << "bad scalar code " << str << std::endl;
	exit(-1);
}

/**
 * Returns the vector field corresponding
 * to the provided C-string:
 *
 * vel - velocity
 *
 * @param[in] s string code of vector field
 * @return specified vector
 */
gp_data::vector gp_data::read_vector(const std::string &str) {

	if      (str=="vel")   return vector::vel;
	else if (str=="none")  return vector::none;

	std::cerr << "bad vector code " << str << std::endl;
	exit(-1);
}

/**
 * Returns the vector-valued function at the provided
 * location for the specified field
 *
 * @param[in] x horizontal coordinate
 * @param[in] y vertical coordinate
 * @param[out] uu horizontal component
 * @param[out] vv vertical component
 */
void gp_data::qp_val(const vector &f,double x,double y,double &uu,double &vv) {

	// grab complex value according to this field
	vcomplex v;
	switch (f) {
		case vector::vel:
			v = br->get_u(x,y);
			if (zero_inf) v -= br->V_inf;
			break;
	}

	dcomplex cu(v.x), cv(v.y);

	// grab scalar value according to this par
	switch (part) {

		case complex_part::real: 
			uu = std::real(cu);
			vv = std::real(cv);
			break;

		case complex_part::imag:
			uu = std::imag(cu);
			vv = std::imag(cv);
			break;
	}
}

/**
 * get the scalar value of the complex field to be plotted
 * as part of a heatmap
 *
 * @param[in] x horiz coordinate
 * @param[in] y vert coordinate
			// get fields
			br->get_u_derivs(x,y,ux_v,ux_x,ux_y,ux_xx,ux_xy,ux_yy);
			br->get_v_derivs(x,y,uy_v,uy_x,uy_y,uy_xx,uy_xy,uy_yy);
			brink_derivs(x,y,ux_v,ux_x,ux_y,ux_xx,ux_xy,ux_yy,
					uy_v,uy_x,uy_y,uy_xx,uy_xy,uy_yy);
 */
double gp_data::hm_val(const scalar &s,double x,double y) {

	// grab complex value according to this field
	dcomplex v;

	// vel values
	dcomplex ux_v,ux_x,ux_y,ux_xx,ux_xy,ux_yy;
	dcomplex uy_v,uy_x,uy_y,uy_xx,uy_xy,uy_yy;
	dcomplex ux_k_v,ux_k_x,ux_k_y,ux_k_xx,ux_k_xy,ux_k_yy;
	dcomplex uy_k_v,uy_k_x,uy_k_y,uy_k_xx,uy_k_xy,uy_k_yy;

	switch (s) {

		// velocity parts
		case scalar::vel_x:

			v = br->get_ux(x,y);
			if (zero_inf) v -= br->V_inf.x;
			break;

		case scalar::vel_y:

			v = br->get_uy(x,y);
			if (zero_inf) v -= br->V_inf.y;
			break;

		// total velocity magnitude
		case scalar::vel_mag:
			{
				double ux(hm_val(scalar::vel_x,x,y));
				double uy(hm_val(scalar::vel_y,x,y));
				return sqrt(ux*ux+uy*uy);
			} break;

		// flow parts
		case scalar::pres:  v = br->get_p(x,y);  break;
		case scalar::vort:  v = br->get_w(x,y);  break;
		case scalar::div:   v = br->get_d(x,y);  break;

		default:
			fputs("bad scalar value\n",stderr);
			exit(-1);
	}

	// grab scalar value according to the part
	switch (part) {
		case complex_part::real: return std::real(v);
		case complex_part::imag: return std::imag(v);
		case complex_part::norm: return std::norm(v);
		case complex_part::arg:  return  std::arg(v);
		default:   fputs("bad part value\n",stderr); exit(-1);
	}
}

/**
 * calculates an array of values representing the arrows
 * of a quiver plot
 *
 * @param[out] spec gnuplot spec string
 * @return the (m*n x 4) array of quiver arrows
 */
mat gp_data::quiver_mat(const vector &f,int m,int n,double sp) {

	// grid spacing 
	double dx((bx-ax)/m), dy((by-ay)/n);

	// data storage
	mat arrows(m*n,4);

	// loop through locations in grid
	double xx,yy,uu,vv,X,Y;
	for(int j=0;j<n;j++) {

		// grab y-coordinate
		yy = ay+(j+0.5)*dy;

		// loop through x-locations
		for(int i=0;i<m;i++) {

			// grab x-coordinate
			xx = ax+(i+0.5)*dx;

			// get velocity vals
			X = rot?yy:xx, Y = rot?xx:yy;
			double &U(rot?vv:uu),&V(rot?uu:vv);
			qp_val(f,X,Y,U,V);

			// velocity magnitude
			double umag = sqrt(uu*uu + vv*vv);
			uu *= sp*dx/umag;
			vv *= sp*dx/umag;

			// get start of vector
			double xs(xx-0.5*uu),ys(yy-0.5*vv);

			// record row
			arrows.row(i+j*m) << xs,ys,uu,vv;
		}
	}

	return arrows;
}

/**
 * calculates an array of values representing the values
 * of a heatmap
 *
 * @param[in] f which scalar field to plot
 * @param[in] m number of cols
 * @param[in] n number of rows
 * @return the (m x n) array of field values
 */
mat gp_data::heatmap_mat(const scalar &f,int m,int n) {

	// grid spacings
	double dx((bx-ax)/m),dy((by-ay)/n);

	// data to output
	mat map(m*n,3);

	// loop through data
	double xx,yy,vv;
	for (int j=0;j<n;j++) {

		// y-coordinate
		yy = ay + (j+0.5)*dy;

		for (int i=0;i<m;i++) {

			// x-coordinate
			xx = ax + (i+0.5)*dx;

			// value
			vv = hm_val(f,rot?yy:xx,rot?xx:yy);

			// store
			map.row(i+j*m) << xx,yy,vv;
		}
	}

	return map;
}
