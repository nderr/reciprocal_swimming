///////////////////////////////////////////////////////////////////////////////
// gp_data.hh                                                                //
//                                                                           //
// code for processing / outputing gnuplot graphics from a brinkman FEM      //
// solution as saved in a binary file                                        //
//                                                                           //
// Nick Derr, 01/31/21                                                       //
///////////////////////////////////////////////////////////////////////////////

#include "brinkman.hh"
#include <eigen3/Eigen/Dense>
#include "gnuplot-iostream.h"

typedef Eigen::ArrayXXf mat;

using namespace gnuplotio;

/**
 * class for processing and outputing gnuplot graphics from
 * a binary save file of the brinkman class
 */
class gp_data {

	public:

	/** scalar parts of complex fields */
	enum class complex_part {real,imag,norm,arg};
	/** possible scalar fields */
	enum class scalar {vel_x,vel_y,vel_mag,pres,vort,div,none};
	/** possible vector fields */
	enum class vector {vel,none};

	/** whether to rotate the grid 90 degrees */
	bool rot;
	/** whether to shift velocities to vanish at the far point */
	bool zero_inf;
	/** horiz grid boundaries */
	double ax,bx;
	/** vert grid boundaries */
	double ay,by;

	/** part of complex field to plot */
	complex_part part;
	/** pointer to two_sphere parameters */
	two_spheres *ts;
	/** pointer to brinkman solution */
	brinkman *br;
	/** pointer to Gnuplot instance */
	Gnuplot *gp;

	// constructors
	gp_data(char *f_in,double ax_,double ay_,double Lx,double Ly,std::string p="r",
		bool zi=true,bool r=true): rot(r),zero_inf(zi),
		ax(ax_),bx(ax+Lx),ay(ay_),by(ay+Ly),part(read_part(p)),
		ts(two_spheres::load_binary(f_in)),br(brinkman::load_binary(f_in)) {}
	~gp_data() {}

	/** helper for memory freeing */
	void free() {delete ts; delete br;}

	// string-to-enum conversions
	static complex_part read_part(const std::string &s);
	static scalar       read_scalar(const std::string &s);
	static vector       read_vector(const std::string &s);

	// field / complex part updates
	void set_part(const complex_part &p) {part=p;}

	// value calculations
	void qp_val(const vector &f,double x,double y,double &u,double &v);
	double hm_val(const scalar &f,double x,double y);

	// data calculations
	mat quiver_mat(const vector &f,int m,int n,double sp=0.8);
	mat heatmap_mat(const scalar &f,int m,int n);

	/**
	 * add the sphere outlines to a plot
	 */
	void add_circles(PlotGroup &plots) {

		// get circle information
		double XL(rot?0:ts->xl),YL(rot?ts->xl:0);
		double XR(rot?0:ts->xr),YR(rot?ts->xr:0);
		mat circs(2,3);
		circs << XL, YL, ts->rl,
	             XR, YR, ts->rr;

		plots.add_plot1d(circs,"u 1:2:(0):3 w circ fc 'black' lc 'black' lw 3");
	}

	/**
	 * add a quiver plot corresponding to the provided vector quantity
	 */
	void add_quiver(vector &f,int mm,int nn,PlotGroup &plots,double sp=0.8) {
		mat dat = quiver_mat(f,mm,nn,sp);
		plots.add_plot1d(quiver_mat(f,mm,nn),"u 1:2:(0):3:4:(0) with vectors arrowstyle 1");
	}

	/**
	 * add a set of streamlines corresponding to the provided vector quantity,
 	 * originating on a m x n grid of points
	 *
	 * @param f which field to plot
	 * @param T duration of time to integrate streamlines
	 * @param tol tolerance for streamline calculations
	 * @param m rows of streamline origination points
	 * @param n columns of streamline origination points
	 * @param plots the group to add plots to
	 * @param cmds commands to have gnuplot execute
	 */
	void add_streamlines(const vector &f,double T,double tol,int m,int n,PlotGroup &plots,std::vector<std::string> &cmds) {

		// set up a set of matrices with streamline data
		std::vector<mat> lines(m*n);
		for (int i=0;i<lines.size();i++) {
			lines[i].resize(10,2);
		}

		// calculate the streamlines
		calculate_streamlines(f,T,tol,m,n,lines);

		// add each streamline
		for (int i=0;i<m*n;i++) {

			// add the plot
			printf("adding plot %d (%d,%d), %d segments\n",i,i%m,i/m,lines[i].rows());
			plots.add_plot1d(lines[i],"u 1:2:(0) w l lw 1 lc rgbcolor 'black'");

			// for each streamline, add an arrowhead
			double len = 0.2;
			if (lines[i].rows() > 1) {

				mat arr;
				arr.resize(1,4);

				int ind = 0;
				double x1 = lines[i](ind  ,0),   y1 = lines[i](ind  ,1);
				double x2 = lines[i](ind+1,0),   y2 = lines[i](ind+1,1);

				double tx(x2-x1), ty(y2-y1);
				double tmag = sqrt(tx*tx+ty*ty);

				tx /= tmag; ty /= tmag;

				double ang = 15;
				double slen = len*cos(ang*180/M_PI);

				// arrow points
				arr(0,0) = x1 - slen*tx;
				arr(0,1) = y1 - slen*ty;
				arr(0,2) = -slen*tx;
				arr(0,3) = -slen*ty;

				std::stringstream vec_plot;
				vec_plot << "u 1:2:(0):3:4:(0) w vec lw 1 lc rgbcolor 'black' head filled size first " << len << "," << ang << ",90 fixed";

				plots.add_plot1d(arr,vec_plot.str());
			}
		}
	}

	/**
	 * Calculate a set of streamlines over an m x n grid
	 *
	 * @param f which field to plot
	 * @param T duration to integrate streamline
	 * @param tol local error for adaptive timestepping
	 * @param m rows of streamline origination points
	 * @param n cols of streamline origination points
	 * @param lines vector of streamlines
	 */
	void calculate_streamlines(const vector &f,double T,double tol,
			int m,int n, std::vector<mat> &lines) {

		// grid spacings
		double dx((bx-ax)/m),dy((by-ay)/n);

		// loop through data
		double xx,yy;

		// RK info
		static const int ns = 4;
		double uu[ns],vv[ns];
		double xi[ns],yi[ns];

		// Butcher tableau
		static const double a[ns][ns] = {
			{    0,    0,    0,    0},
			{ 1./2,    0,    0,    0},
			{    0, 3./4,    0,    0},
			{ 2./9, 1./3, 4./9,    0}
		};
		static const double b3[ns] = {   2./9, 1./3, 4./9,    0};
		static const double b2[ns] = {  7./24, 1./4, 1./3, 1./8};

		// for each streamline
		for (int j=0;j<n;j++) {
			for (int i=0;i<m;i++) {

				printf("(%d/%d) (%d/%d)...",i+1,m,j+1,n);
				fflush(stdout);

				// y-coordinate
				yy = ay + (j+0.5)*dy;

				// x-coordinate
				xx = ax + (i+0.5)*dx;

				double x0(xx),y0(yy);

				// minimum segment length
				double ds = 0.01;
				double t = 0;
				double dt = 0.1;

				// count points
				int pt = 0;
				lines[i+j*m].row(pt++) << xx,yy; 

				bool out=false;

				// that's starting point. now integrate
				while (t < T) {

					double ee = tol+1;
					do {

						// steps of RK
						for (int s=0;s<ns;s++) {

							if (s>0) for (int sp=0;sp<s;sp++) {
								xi[s] = xx - dt*a[s][sp]*uu[sp];
								yi[s] = yy - dt*a[s][sp]*vv[sp];
							} else {
								xi[0] = xx;
								yi[0] = yy;
							}

							double &X(rot?yi[s]:xi[s]), &Y(rot?xi[s]:yi[s]);
							double &U(rot?vv[s]:uu[s]), &V(rot?uu[s]:vv[s]);

							// grab vals for k1
							if (ax < xi[s] && xi[s] < bx && ay < yi[s] && yi[s] < by)
								qp_val(f,X,Y,U,V);
							else
								uu[s] = vv[s] = 0;

							// check nans
							if (std::isnan(uu[s]) || std::isnan(vv[s])) uu[s]=vv[s]=0;
						}


						// grab two estimates for value at next time step
						double x_star(xx),y_star(yy),x_new(xx),y_new(yy);
						for (int s=0;s<ns;s++) {

							x_star -= dt * b2[s] * uu[s];
							y_star -= dt * b2[s] * vv[s];

							x_new  -= dt * b3[s] * uu[s];
							y_new  -= dt * b3[s] * vv[s];
						}

						// compute difference between the two estimates
						// for adaptive timestepping
						ee = sqrt((x_new-x_star)*(x_new-x_star)
								+ (y_new-y_star)*(y_new-y_star));

						fflush(stdout);

						// either accept or reject the new value
						if (ee <= tol) {

							xx = x_new;
							yy = y_new;

							t += dt;
							dt *= 1.1;

						} else {

							dt *= 0.5;

						}

					} while (ee > tol);

				  // store the new value
					double x1 = lines[i+j*m](pt-1,0);
					double y1 = lines[i+j*m](pt-1,1);
					double x2 = xx;
					double y2 = yy;

				  // stop if segment is too short or we reach duration
					double s = sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) );
					if (s > ds || t >= T) {

						if (pt == lines[i+j*m].rows()) {

							mat tmp = lines[i+j*m].eval();
							lines[i+j*m] = mat::Zero(2*pt,2);
							lines[i+j*m].block(0,0,pt,2) = tmp.eval();
						}

						lines[i+j*m].row(pt++) << xx,yy; 
					}

					// check distance from initial point
					double d0 = sqrt((xx-x0)*(xx-x0) + (yy-y0)*(yy-y0));

				  // if we get far from intial point and then get close
				  // again, just stop instead of orbiting
					if (d0 > 3*ds) out = true;
					if (out && d0 < 3*ds) {
						t = T;
					}
				}

				// store, count points, get rid of empty rows
				lines[i+j*m] = lines[i+j*m].block(0,0,pt,2).eval();
				printf("[%d pts]\n",pt);
			}
		}
	}

	/**
	 * Add a heat map of the provided scalar field
	 *
	 * @param f which field to plot
	 * @param mm rows
	 * @param nn columns
	 * @param plots group of plots to add the heatmap to
	 * @param cbr colorbar range
	 */
	void add_heatmap(scalar &f,int mm,int nn,PlotGroup &plots,
		double (&cbr)[2]) {

		// calculate heatmap
		mat dat = heatmap_mat(f,mm,nn);

		// cycle through heatmap, calculating mean
		double mn = 0, norm=0, vl;
		for (int jj=0;jj<nn;jj++) {
			for (int ii=0;ii<mm;ii++) {

				vl = dat(ii+jj*mm,2);

				if (!std::isnan(vl)) {
					mn += vl;
					norm += 1;
				}
			}
		}
		mn /= norm;

		// cycle through heatmap, calculating std dev
		norm = 0;
		double sd = 0;
		for (int jj=0;jj<nn;jj++) {
			for (int ii=0;ii<mm;ii++) {

				vl = dat(ii+jj*mm,2);
				vl *= vl;

				if (!std::isnan(vl)) {
					sd += vl;
					norm += 1;
				}
			}
		}
		sd /= norm;
		sd = sqrt(sd);

		printf("%g +/- %g\n",mn,sd);
		
		// store info in heatmap
		cbr[0] = mn-sd;
		cbr[1] = mn+sd;

		int max_log = -1;
		int min_log = -4;
		
		std::stringstream inst;
		inst << "u 1:2:(0):(abs($3)<10**(" << min_log
			 << ")?0:(sgn($3)*(log10(abs($3))-(" << min_log
			 << ")))) wi imag";
		plots.add_plot1d(dat,inst.str());
	}

	/** add a heatmap to the provided plot group */
	void add_heatmap(scalar &f,int mm,int nn,PlotGroup &plots) {
		double trash[2];
		add_heatmap(f,mm,nn,plots,trash);
	}

	/** start a new gnuplot image */
	void start(bool debug=false) {
		if (debug) gp = new Gnuplot(stdout);
		else       gp = new Gnuplot();
	}

	/** finish a gnuplot image */
	void finish() {delete gp;}

	/** set the plot range */
	void range() {
		*gp << "set xrange [" << ax << ":" << bx << "]" << std::endl;
		*gp << "set yrange [" << ay << ":" << by << "]" << std::endl;
		*gp << "set clabel \'vorticity\'" << std::endl;
	}

	/** pass a command to gnuplot */
	void command(std::string cmd) {*gp << cmd << std::endl;}

	// terminal options
	void term(std::string term,std::string spec="") {
		*gp << "set term " << term << " " << spec << std::endl;
	}
	void term(std::string term,int tm,int tn,std::string spec="") {
		*gp << "set term " << term << " size " << tm
			<< "," << tn << " " << spec << std::endl;
	}

	/**
	 * opens a new pdf for plotting
	 */
	void pdf_open(std::string fname,double m=4.75) {
		char sizes[256];

		double n = m * (by-ay)/(bx-ax)  - 0.8;// - 0.1;  // XXX ND for colorbar

		double xscale=0.075;
		double yscale=xscale * m/n;
		char xscl[256],yscl[256];
		sprintf(xscl,"set lmargin at screen %g",xscale);
		sprintf(yscl,"set bmargin at screen %g",yscale);

		sprintf(sizes," size %gin,%gin",m,n);

		term("pdf","enhanced font ',16' " + std::string(sizes));
		command("set output '" + fname + "'");
		command(xscl);
		command(yscl);
	}
  
	/**
	 * opens a new png for plotting
	 */
	void png_open(std::string fname,int m=550,int n=600) {
		term("pngcairo",m,n,"font ',16'");
		command("set output '" + fname + "'");
	}
	/** close file */
	void close() {command("unset output");}

	/**
	 * create a heatmap with a quiver plot on top
	 */
	void heat_quiver(std::string sca,int sca_m,std::string qui,int qui_m,
			double (&cbrange)[2]) {

		// get from strings
		scalar s = read_scalar(sca);		
		vector v = read_vector(qui);

		// format plots
		range();
		command("set pm3d map");
		command("set pm3d interpolate 4,4");
		command("set style fill solid");
		command("set size ratio -1");
		command("unset key");

		command("set style arrow 1 head filled size 0.12,30,90 lw 1.5 lc 'black'");
		command("set arrow arrowstyle 1");

		// labels
		command("set lmargin at screen 0");
		command("unset xtics");
		command("unset ytics");

		// get row nums
		double ar = (by-ay)/(bx-ax);
		int sca_n = (int) (0.5 + sca_m*ar);
		int qui_n = (int) (0.5 + qui_m*ar);

		// collect plots
		PlotGroup splots = Gnuplot::splotGroup();
		add_heatmap(s,sca_m,sca_n,splots,cbrange);
		add_circles(splots);
		add_quiver(v,qui_m,qui_n,splots);

		// print colorbar range
		printf("%g %g\n",cbrange[0],cbrange[1]);

		// add to gnuplot
		*gp << splots;
	}

	/**
	 * Create a heatmap and overplot with a set of streamlines
	 */
	void heat_stream(std::string sca,int sca_m,std::string qui,int sl_m,
			double sl_T,double sl_tol,
			double (&cbrange)[2]) {

		// get from strings
		scalar s = read_scalar(sca);		
		vector v = read_vector(qui);

		// format plots
		range();
		command("set pm3d map");
		command("set style fill solid");
		command("set size ratio -1");
		command("unset key");

		command("set style arrow 1 head size 0.12,30,90 filled lw 1.5 lc 'black'");
		command("set arrow arrowstyle 1");

		// labels
		command("set lmargin at screen 0");
		command("unset xtics");
		command("unset ytics");

		// get row nums
		double ar = (by-ay)/(bx-ax);
		int sca_n = (int) (0.5 + sca_m*ar);
		int sl_n = (int) (0.5 + sl_m*ar);

		// collect plots
		PlotGroup splots = Gnuplot::splotGroup();
		add_heatmap(s,sca_m,sca_n,splots,cbrange);
		add_circles(splots);

		std::vector<std::string> cmds;
		add_streamlines(v,sl_T,sl_tol,sl_m,sl_n,splots,cmds);

		int min_log = -4;
		int max_log = -1;
		int num_log = 2 * ((max_log - min_log) + 1) + 1;

		// tic labels
		std::stringstream tics;
		tics << "(";
		for(int l = max_log; l >= min_log; l--) {
			tics <<  "\"-10^{" << l << "}\" -(" << (0.5+l-min_log) << "),";
		}
		tics << "0";
		for(int l = min_log; l <= max_log; l++) {
			tics <<  ",\"10^{" << l << "}\"  (" << (0.5+l-min_log) << ")";
		}
		tics << ")";

		// color maps
		std::string low[] = {"\"#6F0695\"", "\"#9630B8\"", "\"#BA58D0\"", 
			"\"#D984DD\"", "\"#F0B3DD\""};
		std::string high[]= {"\"#D6BCE8\"", "\"#AB95F0\"", "\"#7C72E7\"",
			"\"#4B52D0\"", "\"#0F36AB\""};

		double low_v = -(0.5+max_log-min_log);
		double high_v = -low_v;
		double low_z = -0.5;
		double high_z = 0.5;

		*gp << "set palette defined (";
		for (int i=0;i<size(low);i++) {
			*gp << low_v + (low_z-low_v)*i/size(low) << " " << low[i] << ",";
		}

		*gp << low_z << " " << "\"#FFFFFF\"" << "," << high_z
			<< "\"#FFFFFF\"";

		for (int i=0; i<size(high); i++) {
			*gp << "," << high_z + (high_v-high_z)*(i+1)/size(high) << " " << high[i];
		}
		*gp << ")\n";
		*gp << "set cbrange [" << -(0.5+max_log-min_log) << ":" << (0.5+max_log-min_log) << "]; set cbtics " << tics.str() << "\n";
    
		// add to gnuplot
		*gp << splots;

	}
};
