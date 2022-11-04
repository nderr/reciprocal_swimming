#include "quad.hh"

// 1D quadrature

/** constructs the 1D quadrature rule of the requested number of points
 *  using the GNU scientific library implentation, which includes
 *  pre-computed coefficients at low polynomial order */
gauss_quad_1d::gauss_quad_1d(int n_):n(n_),x(new double[n]),w(new double[n]) {

	// get table from GSL Integration
	gsl_integration_glfixed_table *t =
		gsl_integration_glfixed_table_alloc(n);

	// pull out weights/abscissae
	for(int i=0;i<n;i++)
		gsl_integration_glfixed_point(0,1,i,x+i,w+i,t);

	// free table memory
	gsl_integration_glfixed_table_free(t);
}

/** deconstructor */
gauss_quad_1d::~gauss_quad_1d() {
	delete[] w;
	delete[] x;
}

// 2D abstract class

/** constructs a 2D quadrature rule with the desired number of
 *  total points */
gauss_quad_2d::gauss_quad_2d(int N_) : N(N_),x(new double[N]),y(new double[N]),
	w(new double[N]) {}

/** constructs a 2D quadrature rule with the desired number of
 *  total points. only the requested amount of memory is
 *  allocated for the abscissae */
gauss_quad_2d::gauss_quad_2d(int N_,int nx,int ny) : N(N_),x(new double[nx]),
	y(new double[ny]),w(new double[N]) {}

gauss_quad_2d::~gauss_quad_2d() {
	delete[] w;
	delete[] y;
	delete[] x;
}

/** constructs 2D quadrature rule as a tensor product of 1D rule */
void gauss_quad_2d_square_tens::setup_rule() {

	// get 1D Gauss-Legendre rule
	gauss_quad_1d g(n);

	// get coords
	memcpy(x, g.x, n*sizeof(double));
	memcpy(y, x,   n*sizeof(double));

	// tensor product for weights (symmetric)
	for(int j=0;j<n;j++) for(int i=0;i<=j;i++) {

		w[i+n*j] = g.w[i]*g.w[j];
		if (i<j) w[j+n*i] = w[i+n*j];
	}
}

/** constructs the quadrature rule by mapping the abscissae and
 *  weights from the square [0,1]x[0,1] to the standard triangle */
void gauss_quad_2d_tri_tens::setup_rule() {

	// get 2D square Gauss-Legendre (tensor product)
	gauss_quad_2d_square_tens g(n);

	// map such that [x,y] -> [x(1-y),y]

	// y-values stay the same
	memcpy(y, g.y, n*sizeof(double));

	// adjust x-values and weights
	for(int j=0;j<n;j++) for(int i=0;i<n;i++) {
		x[i+j*n] = g.x[i]*(1-g.y[j]);
		w[i+j*n] = g.w[i+j*n]*(1-g.y[j]);
	}
}
