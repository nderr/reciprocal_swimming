#include "brinkman.hh"
#include <gsl/gsl_sf_legendre.h>
#include <cmath>

int main(int argc,char **v) {

	// file name
	char *f_in = v[1];

	// sphere index
	int sn = atoi(v[2]);


	// fill in spherical harmonic values
	int n = 1000;
	double dx = 0.01;

	two_spheres *ts = two_spheres::load_binary(f_in);
	brinkman    *br = brinkman::load_binary(f_in);

	double xc = sn == 1 ? ts->xl : ts->xr;
	double r0 = sn == 1 ? ts->rl : ts->rr;
	dcomplex U = sn == 1 ? ts->V_l.x : ts->V_r.x;

	const double eps = 1e-12;
	double rr = r0+eps;

	bool left = sn==1;

	for (int ix=0;ix<n; ix++) {

		double xx_1 = left?(ts->xl-rr-ix*dx):(ts->xr+rr+ix*dx);
		double yy_1 = eps;

		double xx_2 = left?ts->xl:ts->xr;
		double yy_2 = rr + ix*dx;

		// grab solution values
		dcomplex ux,ux_x,ux_y,ux_xx,ux_xy,ux_yy;
		dcomplex uy,uy_x,uy_y,uy_xx,uy_xy,uy_yy;
		dcomplex p,p_x,p_y,p_xx,p_xy,p_yy;

		br->get_u_derivs(xx_1,yy_1,ux,ux_x,ux_y,ux_xx,ux_xy,ux_yy);
		br->get_v_derivs(xx_2,yy_2,uy,uy_x,uy_y,uy_xx,uy_xy,uy_yy);
		br->get_p_derivs(xx_2,yy_2,p,p_x,p_y,p_xx,p_xy,p_yy);


		//printf("%g %g\n",tt, std::real(I * 0.5 * std::conj(U) *tr));
		printf("%g %g %g %g %g %g %g\n", rr+ix*dx,
				std::real(ux), std::imag(ux),
				std::real(uy), std::imag(uy),
				std::real(p),  std::imag(p));


		//f_out += std::real(I * 0.5 * std::conj(U) * tr);
	}
}
