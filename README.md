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
