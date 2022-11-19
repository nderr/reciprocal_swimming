#include <cstdio>
#include <cstdlib>
#include "petscksp.h"
#include "arg_manager.hh"
#include "iostream"

#include "brinkman.hh"

static char help[] = "testing Petsc for solving Stokes";

// forward declare arg proc
void process_args(arg_manager &args);

int main(int argc,char **argv) {

	// MPI and Initialization
	PetscInitialize(&argc,&argv,0,help);

	// set up arg input
	arg_manager args(argc,argv,2);
	process_args(args);

	// grab args
	int m = args.get_int(0);
	int p = args.get_int(1);

	// inertia
	double M_sq  = args.get_double('M');
	double alpha = sqrt(M_sq * args.get_double('f'));

	// geometry
	double R     = args.get_double('R');
	double r     = args.get_double('r');
	double d     = args.get_double('s');

	// viscosity
	double mu    = args.get_double('m');

	// densities
	double rho_l = M_sq * args.get_double('e');
	double rho_r = M_sq * args.get_double('d');

	// instantiate and set problem info
	two_spheres sim_brink(p,m,m/2,d,R,r,mu,alpha,-1,1);
	sim_brink.sol_type = petsc_solve::method::lu;
	sim_brink.set_densities(rho_l,rho_r);

	// buffer for strings
	char fn[256];

	// solve subject to equal internal forces
	sim_brink.solve_force_balance();

	// normalize and re-solve
	sim_brink.norm_Vx();
	sim_brink.solve();

	// output binary file
	sprintf(fn,"%s_b.bin",args.get_string('o'));
	sim_brink.save_binary(fn);

	// steady stokes solution
	two_spheres sim_stokes(&sim_brink,0,args.get_int('p'));
	sim_stokes.sol_type = petsc_solve::method::lu;

	// stokes
	sim_stokes.solve_force_free();

	// output binary file
	sprintf(fn,"%s_s.bin",args.get_string('o'));
	sim_stokes.save_binary(fn);

	// grab swim speed
	dcomplex Us = -sim_stokes.V_inf.x;

	// print output
	mesg("%010.8f %010.8f %010.8f %010.8f\n",M_sq,rho_l,rho_r,std::real(Us));

	// clean up PETSC datatypes
	PetscFinalize();
};

/**
 * Processes command-line args with the arg_manager class
 */
void process_args(arg_manager &args) {

	// grab element grid / poly order
	args.add_int("m","horiz grid");
	args.add_int("p","poly order");

	// opts for geometry
	args.add_double('s',3,"separation","distance between sphere centers");
	args.add_double('R',1,"left r","radius of left sphere");
	args.add_double('r',0.5,"right r","radius of right sphere");

	// inertia and density
	args.add_double('M',1,"M_sqrd","Riley's M^2 parameter");
	args.add_double('f',1,"fluid inertia","M^2 parameter of fluid in units of "
			"left sphere density");
	args.add_double('m',1,"mu","viscosity");
	args.add_double('d',1,"right density","right density (in units of left density)");
	args.add_double('e',1,"left density","left density (in units of left density)");

	// problem types
	args.add_int('p',0,"part of steady sol","which part of stead sol\n"
			"0: total soln\n"
			"1: only bvel\n"
			"2: only bulk force");

	// file output
	args.add_string('o',"out","outfile","binary output filename");

	args.process_args();
}
