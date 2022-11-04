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

	// grab args
	arg_manager args(c,v,1);

	gp_data gp = process_args(args);

	gp.start();
	//gp.png_open(std::string(args.get_string('o'))+".png");
	gp.pdf_open(std::string(args.get_string('o'))+".pdf");

//	gp.term("six",900,1100);

	// which fields to plot
	std::string heat = args.get_string('s');
	std::string quiv = args.get_string('q');

	// get resolutions
	int heat_m = args.get_int('P');
	int quiv_m = args.get_int('A');

	double cbrange[2];

	std::string cbr = args.get_string('c');
	gp.command("set cbrange "+cbr);
	std::string ttl = args.get_string('t');
	gp.command("set title '" + ttl + "'");

	double sl_T  = args.get_double('T');
	double sl_dt = args.get_double('d');


	//gp.command("unset colorbox");
	//gp.command("set cblabel 'vorticity' offset -3");
	gp.heat_stream(heat,heat_m,quiv,quiv_m,
			sl_T,sl_dt,cbrange);


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

	// which field to output
	args.add_string('s',"none","heatmap field",
			"scalar simulation field to output");
	args.add_string('q',"none","quiver field",
			"vector field to output as quiver plot");
	args.add_string('l',"none","streamline field",
			"vector field to output as streamlines");
	args.add_string('p',"r",   "complex part",
			"scalar part of complex field");

	args.add_string('o',"out","outfile","file to write to (ext appended)");

	// check for output
	

	// heatmap resolution
	args.add_int('P',400,"x-pixels","number of pixels in x-direction for heat map");
	args.add_int('A',36,"x-arrows","number of arrows in x-direction for quiver plot");

	// grid extent (w/in brinkman binary)
	args.add_double('X',9,"horiz extent","horizontal domain plotting extent");
	args.add_double('Y',12,"vert extent","vertical domain plotting extent");
	args.add_double('x',-0.5*args.get_double('X'),"left edge","left edge of domain");
	args.add_double('y',-0.5*args.get_double('Y'),"bottom edge","bottom edge of domain");
	args.add_double('T',1,"streamline time","duration of streamline integration");
	args.add_double('d',0.1,"streamline time disc","time discretation of streamline integration");

	args.add_switch('I',"nonzero V_inf","don't subtract off the velocity at the far point");
	args.add_switch('R',"rotate","whether to rotate 90 degrees");

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
