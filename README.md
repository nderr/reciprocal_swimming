This repository contains the code used to generate the data and figures presented in the *Journal of Fluid Mechanics* publication [Reciprocal swimming at intermediate Reynolds number](https://arxiv.org/abs/2202.03669).

### Dependencies

Dependencies for simulation
- [PETSc](https://petsc.org/release/) with [MPI](https://www.open-mpi.org/) for parallelization
- [GNU Scientific Library](https://www.gnu.org/software/gsl/)
- [BLAS](https://www.openblas.net/)
- C++-17

Dependencies for visualization
- [gnuplot-iostream](https://github.com/dstahlke/gnuplot-iostream)
- [Boost](https://www.boost.org/)
  - [Boost.Iostreams](https://www.boost.org/doc/libs/1_80_0/libs/iostreams/doc/index.html)
  - [Boost.Filesystem](https://www.boost.org/doc/libs/1_78_0/libs/filesystem/doc/index.htm)

### Configuration
1. Clone the repository
2. Copy the file config.mk.in to config.mk
3. Edit the contents of config.mk to correspond to the locations on your system of the above dependencies

### Building
After navigating into the repository and creating/editing config.mk as described above, "make" will build the executables. If you're only interested in building the simulation code, you can run "make solves" to build the solver without the visualization tools (avoiding the need to install their dependencies.) Similarly, "make stats" will build only the visualization tools.

### Running
To solve a system with separation distance S, left sphere radius lr, right sphere radius rr, and inertial parameter M^2, run

mpirun -np num_procs swim_steady -S S -R lr -r rr -M M^2 -o output_directory m p

to solve the system over a 2m x m bispherical grid of finite elements. The velocity elements are of polynomial order p, and the pressure elements are of order p-1. The code will print to stdout:
- the provided M^2 valu
- the density of the left sphere
- the density of the right sphere
- the resulting steady swim speed

Run "./swim_steady -h" for a short help message and "./swim_steady -H" for more documentation on option parameters. Within the output directory, the Brinkman and Stokes flow fields will also be saved as binary files. The executable "./gen_gp" can be used to read in these files and construct heat maps of the Stokes vorticity field with overlaid streamlines.

### Contact
For questions about the code or the 
