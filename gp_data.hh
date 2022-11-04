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
	enum class scalar {vel_x,vel_y,vel_mag,pres,vort,div,none,err,err_x,err_y,err_xx,err_xy,err_yy};
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

	gp_data(char *f_in,double ax_,double ay_,double Lx,double Ly,std::string p="r",
		bool zi=true,bool r=true): rot(r),zero_inf(zi),
		ax(ax_),bx(ax+Lx),ay(ay_),by(ay+Ly),part(read_part(p)),
		ts(two_spheres::load_binary(f_in)),br(brinkman::load_binary(f_in)) {}
	~gp_data() {}

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

	void add_circles(PlotGroup &plots) {

		// get circle information
		double XL(rot?0:ts->xl),YL(rot?ts->xl:0);
		double XR(rot?0:ts->xr),YR(rot?ts->xr:0);
		mat circs(2,3);
		circs << XL, YL, ts->rl,
	             XR, YR, ts->rr;

		plots.add_plot1d(circs,"u 1:2:(0):3 w circ fc 'black' lc 'black' lw 3");
	}

	void add_quiver(vector &f,int mm,int nn,PlotGroup &plots,double sp=0.8) {
		mat dat = quiver_mat(f,mm,nn,sp);
		plots.add_plot1d(quiver_mat(f,mm,nn),"u 1:2:(0):3:4:(0) with vectors arrowstyle 1");
	}

	void add_streamlines(const vector &f,double T,double tol,int m,int n,PlotGroup &plots,std::vector<std::string> &cmds) {

		/*
		double df = T/nf;

		int per;
		if (df < dt) {
			dt = df;
			per = 1;
		} else {

			double per_d = df/dt;
			per = int(per_d);
			double rem = df - per*dt;

			dt -= rem/(++per);

		}
		*/


		std::vector<mat> lines(m*n);
		for (int i=0;i<lines.size();i++) {
			lines[i].resize(10,2);
		}

		calculate_streamlines(f,T,tol,m,n,lines);

		for (int i=0;i<m*n;i++) {

			printf("adding plot %d (%d,%d), %d segments\n",i,i%m,i/m,lines[i].rows());

			//printf("%d %d\n",lines[i].rows(),lines[i].cols());
//			printf("%g %g\n%g %g\n...",lines[i](0,0),lines[i](0,1),lines[i](1,0),lines[i](1,1));
			plots.add_plot1d(lines[i],"u 1:2:(0) w l lw 1 lc rgbcolor 'black'");

			double len = 0.2;

			if (lines[i].rows() > 1) {

				mat arr;
				arr.resize(1,4);


				int ind = 0;// lines[i].rows()-2; //(lines[i].rows()-1) / 2;
				double x1 = lines[i](ind  ,0),   y1 = lines[i](ind  ,1);
				double x2 = lines[i](ind+1,0),   y2 = lines[i](ind+1,1);

				double tx(x2-x1), ty(y2-y1);
				double tmag = sqrt(tx*tx+ty*ty);

				tx /= tmag; ty /= tmag;

				double ang = 15;
				double slen = len*cos(ang*180/M_PI);


				arr(0,0) = x1 - slen*tx;
				arr(0,1) = y1 - slen*ty;
				arr(0,2) = -slen*tx;
				arr(0,3) = -slen*ty;

				std::stringstream vec_plot;
				vec_plot << "u 1:2:(0):3:4:(0) w vec lw 1 lc rgbcolor 'black' head filled size first " << len << "," << ang << ",90 fixed";

				plots.add_plot1d(arr,vec_plot.str());

				//printf("%g %g %g %g\n",x1,y1,x2-x1,y2-y1);

			}

			//cmds.emplace_back(arr_cmd.str());

		}


	}

	void calculate_streamlines(const vector &f,double T,double tol,
			int m,int n, std::vector<mat> &lines) {

		// grid spacings
		double dx((bx-ax)/m),dy((by-ay)/n);

		// loop through data
		double xx,yy;

		static const int ns = 4;

		double uu[ns],vv[ns];
		double xi[ns],yi[ns];

		static const double a[ns][ns] = {
			{    0,    0,    0,    0},
			{ 1./2,    0,    0,    0},
			{    0, 3./4,    0,    0},
			{ 2./9, 1./3, 4./9,    0}
		};
		static const double b3[ns] = {   2./9, 1./3, 4./9,    0};
		static const double b2[ns] = {  7./24, 1./4, 1./3, 1./8};

//		static const double a[ns][ns] = {
//			{    0,    0,    0},
//			{ 1./3,    0,    0},
//			{    0, 2./3,    0}
//		};
//		static const double b3[ns] = {1./4,   0, 3./4}, b2[ns] = {};

		for (int j=0;j<n;j++) {

			for (int i=0;i<m;i++) {

				printf("(%d/%d) (%d/%d)...",i+1,m,j+1,n);
				fflush(stdout);

				// lines[i+j*m].resize(nf,2);

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


						double x_star(xx),y_star(yy),x_new(xx),y_new(yy);
						for (int s=0;s<ns;s++) {

							x_star -= dt * b2[s] * uu[s];
							y_star -= dt * b2[s] * vv[s];

							x_new  -= dt * b3[s] * uu[s];
							y_new  -= dt * b3[s] * vv[s];
						}

						ee = sqrt((x_new-x_star)*(x_new-x_star)
								+ (y_new-y_star)*(y_new-y_star));

						//printf("dt=%g: (%.8g,%.8g) - (%.8g,%.8g) - (%.8g,%.8g)\n",
						//			dt,xx,yy,x_star,y_star,x_new,y_new);

						//printf("[%g->%g]...",dt,ee);
						fflush(stdout);

						if (ee <= tol) {

							xx = x_new;
							yy = y_new;

							t += dt;
							dt *= 1.1;

						} else {

							dt *= 0.5;

						}

					} while (ee > tol);

					double x1 = lines[i+j*m](pt-1,0);
					double y1 = lines[i+j*m](pt-1,1);
					double x2 = xx;
					double y2 = yy;

					//printf("dt=%g: %g %g %g %g\n",dt,x1,y1,x2,y2);

					double s = sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) );

					//printf("%g\n",s);

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

					if (d0 > 3*ds) out = true;

					if (out && d0 < 3*ds) {
						t = T;
					}
				}

				lines[i+j*m] = lines[i+j*m].block(0,0,pt,2).eval();

				printf("[%d pts]\n",pt);
			}
		}
	}

	void add_heatmap(scalar &f,int mm,int nn,PlotGroup &plots,
		double (&cbr)[2]) {

		mat dat = heatmap_mat(f,mm,nn);

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

	void add_heatmap(scalar &f,int mm,int nn,PlotGroup &plots) {
		double trash[2];
		add_heatmap(f,mm,nn,plots,trash);
	}

	void start(bool debug=false) {
		if (debug) gp = new Gnuplot(stdout);
		else       gp = new Gnuplot();
	}

	void finish() {delete gp;}

	void range() {
		*gp << "set xrange [" << ax << ":" << bx << "]" << std::endl;
		*gp << "set yrange [" << ay << ":" << by << "]" << std::endl;
		*gp << "set clabel \'vorticity\'" << std::endl;
	}
	void command(std::string cmd) {*gp << cmd << std::endl;}
	void term(std::string term,std::string spec="") {
		*gp << "set term " << term << " " << spec << std::endl;
	}
	void term(std::string term,int tm,int tn,std::string spec="") {
		*gp << "set term " << term << " size " << tm
			<< "," << tn << " " << spec << std::endl;
	}

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
		//command("set palette defined (0. \"#6F0695\", 1./10. \"#9630B8\",2./10. \"#BA58D0\", 3./10. \"#D984DD\",4./10. \"#F0B3DD\", 5./10. \"#FFFFFF\" , 6./10. \"#D6BCE8\" , 7./10. \"#AB95F0\",8./10. \"#7C72E7\",9./10. \"#4B52D0\",1. \"#0F36AB\")");
//		command("set tmargin at screen 0.95");
//		command("set rmargin at screen 0.8");
//		command("set lmargin at screen 0.05");
	}
	void png_open(std::string fname,int m=550,int n=600) {
		term("pngcairo",m,n,"font ',16'");
		command("set output '" + fname + "'");
	}
	void close() {command("unset output");}

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
//		command("set style arrow 1 nohead lc 'white'");
		command("set arrow arrowstyle 1");

		// labels
//		command("set xlabel 'x'");
//		command("set ylabel 'y'");
//		command("set title 'velocity magnitude'");
//		command("set cblabel 'velocity magnitude'");
		command("set lmargin at screen 0");
//		command("set rmargin at screen 1");
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

		printf("%g %g\n",cbrange[0],cbrange[1]);




		// add to gnuplot
		*gp << splots;
	}

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
//		command("set style arrow 1 nohead lc 'white'");
		command("set arrow arrowstyle 1");

		// labels
//		command("set xlabel 'x'");
//		command("set ylabel 'y'");
//		command("set title 'velocity magnitude'");
//		command("set cblabel 'velocity magnitude'");
		command("set lmargin at screen 0");
//		command("set rmargin at screen 1");
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
		//add_quiver(v,qui_m,qui_n,splots);

		std::vector<std::string> cmds;
		add_streamlines(v,sl_T,sl_tol,sl_m,sl_n,splots,cmds);

		//*gp << "set cbrange [" << cbrange[0]/10 <<
		//	":" << cbrange[1]/10 << "]; unset cblabel\n";
		//*gp << "set cbrange [*:*]\n";
		//
		int min_log = -4;
		int max_log = -1;
		int num_log = 2 * ((max_log - min_log) + 1) + 1;

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


		//std::string cb_tics = "(\"-10^{-1}\" -4,\"-10^{-3}\" -2,0,\"10^{-3}\" 2,\"10^{-1}\" 4)";
		//

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

		//*gp << "set palette defined (0. \"#6F0695\", 1./10. \"#9630B8\",2./10. \"#BA58D0\", 3./10. \"#D984DD\",4./10. \"#F0B3DD\", 5./10. \"#FFFFFF\" , 6./10. \"#D6BCE8\" , 7./10. \"#AB95F0\",8./10. \"#7C72E7\",9./10. \"#4B52D0\",1. \"#0F36AB\")");
		*gp << "set cbrange [" << -(0.5+max_log-min_log) << ":" << (0.5+max_log-min_log) << "]; set cbtics " << tics.str() << "\n";
		// add to gnuplot
		*gp << splots;

		//for (std::string c : cmds) {
		//	*gp << c;
		//}
	}
};
