///////////////////////////////////////////////////////////////////////////////
// gen_gp.cc                                                                 //
//                                                                           //
// generates gnuplot input files from a brinkman solution binary             //
//                                                                           //
// Nick Derr, 1/31/21                                                        //
// cleanup 1/29/21                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "arg_manager.hh"
#include "gp_data.hh"

gp_data process_args(arg_manager &args);

int main(int c,char **v) {

	PetscInitialize(&c,&v,0,NULL);

	// grab args and read in data
	arg_manager args(c,v,1);
	gp_data gp = process_args(args);

	// open a PDF
	gp.start();
	gp.pdf_open(std::string(args.get_string('o'))+".pdf");

	// grab which fields to plot
	std::string heat = args.get_string('s');
	std::string quiv = args.get_string('q');

	// get resolution and # of streamlines
	int heat_m = args.get_int('P');
	int quiv_m = args.get_int('A');

	// grab color bar range and titles
	double cbrange[2];
	std::string cbr = args.get_string('c');
	gp.command("set cbrange "+cbr);
	std::string ttl = args.get_string('t');
	gp.command("set title '" + ttl + "'");

	// grab duration and dt for streamline
	// integration
	double sl_T  = args.get_double('T');
	double sl_dt = args.get_double('d');

	// make heatmap + streamline
	gp.heat_stream(heat,heat_m,quiv,quiv_m,
			sl_T,sl_dt,cbrange);

	// finish plot
	gp.close();
	gp.finish();
	gp.free();

	PetscFinalize();
}

/**
 * processes command-line arguments
 *
 * @param[out] args command-line argument manager
 */
gp_data process_args(arg_manager &args) {

	// file input / output
	args.add_string("infile", "binary file from PETSC solve to process");

	// which fields to output
	args.add_string('s',"none","heatmap field",
			"scalar simulation field to output");
	args.add_string('q',"none","streamline field",
			"vector field to output as streamline plot");
	args.add_string('p',"r",   "complex part",
			"which scalar part of complex field");

	args.add_string('o',"out","outfile","file to write to (ext appended)");

	// heatmap resolution
	args.add_int('P',400,"x-pixels","number of pixels in x-direction for heat map");
	args.add_int('A',36," x-arrows","number of streamlines in x-direction for streamline plot");

	// grid extent (w/in brinkman binary)
	args.add_double('X',9,"horiz extent","horizontal domain plotting extent");
	args.add_double('Y',12,"vert extent","vertical domain plotting extent");
	args.add_double('x',-0.5*args.get_double('X'),"left edge","left edge of domain");
	args.add_double('y',-0.5*args.get_double('Y'),"bottom edge","bottom edge of domain");

	// streamline integration info
	args.add_double('T',1,"streamline time","duration of streamline integration");
	args.add_double('d',0.1,"streamline time disc","time discretation of streamline integration");

	// rotation and far-field velocities
	args.add_switch('I',"nonzero V_inf","don't subtract off the velocity at the far point");
	args.add_switch('R',"rotate","whether to rotate 90 degrees");

	// colorbar and title
	args.add_string('c',"[*:*]","colorbar range","range for gnuplot \"set cbrange\"");
	args.add_string('t',"","plot title","title of plot");

	args.process_args();

	// in-file
	char *f_in = args.get_string(0);

	// pull in simulation domain range
	const double ax=args.get_double('x');
	const double ay=args.get_double('y');
	const double Lx=args.get_double('X');
	const double Ly=args.get_double('Y');

	// complex part
	std::string cp = args.get_string('p');

	// whether to subtract off V_inf
	bool sub_inf = !args.present('I');
	bool rot = !args.present('R');

	return gp_data(f_in,ax,ay,Lx,Ly,cp,sub_inf,rot);
}
