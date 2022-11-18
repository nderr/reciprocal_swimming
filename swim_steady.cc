#include <cstdio>
#include <cstdlib>
#include "petscksp.h"
#include "omp.h"
#include "arg_manager.hh"
#include "iostream"

#include "brinkman.hh"

static char help[] = "testing Petsc for solving Stokes";

// forward declare arg proc
void process_args(arg_manager &args);

int main(int argc,char **argv) {

	// MPI and Initialization
	PetscInitialize(&argc,&argv,0,help);

	// for now don't do OMP
#if _OPENMP
	omp_set_num_threads(1);
#endif

	// set up arg intput
	arg_manager args(argc,argv,2);
	process_args(args);

	// grab args
	int m = args.get_int(0);
	int p = args.get_int(1);

	double M_sq  = args.get_double('M');

	double alpha = sqrt(M_sq * args.get_double('f'));

	double R     = args.get_double('R');
	double r     = args.get_double('r');
	double d     = args.get_double('s');

	double mu    = args.get_double('m');

	double rho_l = M_sq * args.get_double('e');
	double rho_r = M_sq * args.get_double('d');


	// scale lengths
	double L_sca = R > r ? R : r;
	R /= L_sca; r /= L_sca; d /= L_sca;

	// instantiate problem info and stokes class
	two_spheres sim_brink(p,m,m/2,d,R,r,mu,alpha,-1,1);

	int levels = 1, mm=m;
	while (mm >4) {levels++; mm>>=1;}

	switch (args.get_int('t')) {
		case 0: sim_brink.sol_type = petsc_solve::method::lu;   break;
		case 1: sim_brink.sol_type = petsc_solve::method::mg;   break;
		case 2: sim_brink.sol_type = petsc_solve::method::schur;break;
	}

	sim_brink.set_densities(rho_l,rho_r);

	// 


//	sim_brink.solve_petsc_inplace();
//	sim_brink.solve_petsc_dmda();

	char fn[256];
//	sprintf(fn,"%s_b.bin",args.get_string('o'));
//	sim_brink.save_binary(fn);
//	PetscFinalize();
//	exit(0);

//	if (!args.present('D')) sim_brink.solve_petsc_dmda();
//	else                    sim_brink.solve_petsc_ksp();

	//*
	if (args.present('E')) {
		sim_brink.solve();
		sprintf(fn,"%s_b.bin",args.get_string('o'));
		sim_brink.save_binary(fn);
		sim_brink.ux_r_quad("quad_x");
		sim_brink.ux_r_line("line_x");
		bail(0,"done for now\n");
	}
// 	*/




	sim_brink.solve_force_balance();



//	mesg("normalize\n");
	sim_brink.norm_Vx();

	sim_brink.solve();



//	int rank=0;
//	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

//	char fn[256];
	sprintf(fn,"%s_b.bin",args.get_string('o'));
	sim_brink.save_binary(fn);
	sim_brink.ux_r_quad("quad_x");
	sim_brink.ux_r_line("line_x");

//	PetscFinalize();
//	exit(0);

	// steady stokes solution
	two_spheres sim_stokes(&sim_brink,0,args.get_int('p'));
	switch (args.get_int('t')) {
		case 0: sim_stokes.sol_type = petsc_solve::method::lu;   break;
		case 1: sim_stokes.sol_type = petsc_solve::method::mg;   break;
		case 2: sim_stokes.sol_type = petsc_solve::method::schur;break;
	}

	// stokes
	sim_stokes.solve_force_free(args.present('F'));


//	mesg("%g %g\n",std::real(drag_C),std::imag(drag_C));

	sprintf(fn,"%s_s.bin",args.get_string('o'));
	sim_stokes.save_binary(fn);

	dcomplex Ul =  sim_brink.V_l.x;
	dcomplex Ur =  sim_brink.V_r.x;
	dcomplex Us = -sim_stokes.V_inf.x;

	// felderhof
//	double omega = 20*M_PI;
//	double ss = sqrt(R*R*omega*Re / (2*mu));

	dcomplex F_drag_1 = sim_brink.drag_force(true),
			 F_drag_2 = sim_brink.drag_force(false),
			 F_net_1  = sim_brink.net_force(true),
			 F_net_2  = sim_brink.net_force(false);
//			 F_drag = std::real(sim_stokes.drag_force(true)+sim_stokes.drag_force(false));

	dcomplex F_s_l = F_net_1 - F_drag_1;
	dcomplex F_s_r = F_net_2 - F_drag_2;

	dcomplex drag_C = sim_stokes.drag_coeff();
	//mesg("%g %g\n",std::real(drag_C),std::imag(drag_C));





//	double Re = alpha*alpha;

	/*
	mesg("%010.8f %010.8f %010.8f %010.8f %010.8f %010.8f %010.8f %010.8f %010.8f %010.8f\n",
			M_sq,std::real(Us),
			 std::real(F_drag_1),  std::imag(F_drag_1),
			 std::real(F_drag_2),  std::imag(F_drag_2),
             std::real(F_net_1),   std::imag(F_net_1),
			 std::real(F_net_2),   std::imag(F_net_2));
			 */


	mesg("%010.8f %010.8f %010.8f %010.8f\n",
			M_sq,rho_l,rho_r,std::real(Us)
			);

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
//	args.add_double('a',0.2,"alpha","Riley's |M| parameter");
	args.add_double('M',1,"M_sqrd","Riley's M^2 parameter");

	args.add_double('f',1,"fluid inertia","M^2 parameter of fluid in units of "
			"left sphere density");

	args.add_double('m',1,"mu","viscosity");

	args.add_double('d',1,"right density","right density (in units of left density)");
	args.add_double('e',1,"left density","left density (in units of left density)");

	args.add_switch('F',"approx drag coeff","use approximate drag coefficient to match Felderhof");

	args.add_switch('E',"early stop","only do one solve (for debugging e.g.)");
	args.add_switch('D',"don't use dmda","don't use dmda");
	args.add_int('p',0,"part of steady sol","which part of stead sol\n"
			"0: total soln\n"
			"1: only bvel\n"
			"2: only bulk force");


	args.add_int('t',0,"solver type","which solver to use\n0: direct (LU)\n1:"
			" GMRES (ILU PC)\n2: multigrid (GMRES smoothing)\n3: Schur fieldsplit");

	args.add_string('o',"out","outfile","binary output filename");

	args.process_args();
}
