#include "brinkman.hh"
#include <gsl/gsl_sf_legendre.h>
#include <cmath>

int main(int argc,char **v) {

	// file name
	char *f_in = v[1];

	// sphere index
	int sn = atoi(v[2]);


	// fill in spherical harmonic values
	int nt = 10;
	double dt = M_PI / (double) nt;

	two_spheres *ts = two_spheres::load_binary(f_in);
	brinkman    *br = brinkman::load_binary(f_in);

	double xc = sn == 1 ? ts->xl : ts->xr;
	double r0 = sn == 1 ? ts->rl : ts->rr;
	dcomplex U = sn == 1 ? ts->V_l.x : ts->V_r.x;

	const double eps = 1e-12;
	double rr = r0+eps;

	double f_out=0;

	double tot=0;

	for (int it=0;it<nt; it++) {

		// convert polar to cartesian
		double tt = (it+0.5)*dt;
		double xx(rr * cos(tt) + xc), yy(rr * sin(tt));

		// grab solution values
		dcomplex ux,ux_x,ux_y,ux_xx,ux_xy,ux_yy;
		dcomplex uy,uy_x,uy_y,uy_xx,uy_xy,uy_yy;
		dcomplex p,p_x,p_y,p_xx,p_xy,p_yy;

		br->get_u_derivs(xx,yy,ux,ux_x,ux_y,ux_xx,ux_xy,ux_yy);
		br->get_v_derivs(xx,yy,uy,uy_x,uy_y,uy_xx,uy_xy,uy_yy);
		br->get_p_derivs(xx,yy,p,p_x,p_y,p_xx,p_xy,p_yy);

		// grab normal vector and derivs
		double nx = cos(tt);
		double ny = sin(tt);
		double nx_x =  ny*ny/rr;
		double ny_x = -nx*ny/rr;

		dcomplex s_xx   = 2 * ux_x - p;
		dcomplex s_xx_x = 2 * ux_xx - p_x;
		dcomplex s_xy   = ux_y  + uy_x;
		dcomplex s_xy_x = ux_xy + uy_xx;

		dcomplex tr = s_xx_x * nx +
                      s_xy_x * ny +
					  s_xx   * nx_x +
					  s_xy   * ny_x;

		double d_tr = std::real(I * 0.5 * U * std::conj(tr));

		// convert polar to cartesian
		xx = (rr) * cos(tt) + xc;
		yy = (rr) * sin(tt);

		ux = br->get_ux(xx,yy);

		//printf("%g %g\n",tt, std::real(I * 0.5 * std::conj(U) *tr));
		printf("%g %g (%g) %g %g\n",tt, d_tr, d_tr * dt * 2 * M_PI * sin(tt) * rr * rr,
				std::real(ux_x),std::imag(ux_x));

		tot += d_tr * dt * 2 * M_PI * sin(tt) * rr * rr;

		//f_out += std::real(I * 0.5 * std::conj(U) * tr);
	}

	printf("total: %g\n",tot);

}
