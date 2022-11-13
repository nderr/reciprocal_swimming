This repository contains the code used to generate the data and figures presented in the *Journal of Fluid Mechanics* publication "Reciprocal swimming at intermediate Reynolds number."

### Dependencies

Dependencies for simulation
- [PETSc](https://petsc.org/release/)
- [GNU Scientific Library](https://www.gnu.org/software/gsl/)
- BLAS implementation, e.g. [OpenBLAS](https://www.openblas.net/)
- C++-17

Dependencies for visualization
- [gnuplot-iostream](https://github.com/dstahlke/gnuplot-iostream)
- [Boost](https://www.boost.org/)
  - [Boost.Iostreams](https://www.boost.org/doc/libs/1_80_0/libs/iostreams/doc/index.html)
  - [Boost.Filesystem](https://www.boost.org/doc/libs/1_78_0/libs/filesystem/doc/index.htm)

### Configuration
1. Clone the repository
2. Copy the file config.mk.in to config.mk
3. Edit the contents of config.mk to correspond to the the system locations of the above dependencies

### Building
After navigating into the repository and creating/editing config.mk as described above, "make" will build the executables. If you're only interested in building the simulation code, you can type "make solves" to build the solver without the visualization tools (avoiding the need to install their dependencies.) Similarly, "make stats" will build only the visualization tools.

### Running
(More detail on running the code will be added soon.)

To solve a system with separation distance S, left sphere radius lr, right sphere radius rr, and inertial parameter M^2, run

mpirun -np num_procs swim_steady -S S -R lr -r rr -M M^2 -o output_directory m p

This will solve the system over a 2m x m grid of finite elements of polynomial order p for the velocities and p-1 for the pressures.
