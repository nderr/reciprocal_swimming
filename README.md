# Reciprocal swimming at intermediate Reynolds number
This repository contains C++ code used to calculate the stready flow field and swim speed of a simple reciprocal swimmer in the following publication:

N.J. Derr, T. Dombrowski, C.H. Rycroft, and D. Klotsa. [Reciprocal swimming at intermediate Reynolds number](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/reciprocal-swimming-at-intermediate-reynolds-number/DB26C467722B4DE45441ACBB92FD34B3). Journal of Fluid Mechanics 952, A8 (2022). [[arxiv]](https://arxiv.org/abs/2202.03669)

The repository's initial commit is the version of the code base used to generate data and figures on Linux architectures, compiled using GCC 12.1.0. The present version reflects edits to improve documentaion and readability and remove spurious features. 

## Dependencies

### Simulation
- [PETSc](https://petsc.org/release/) with [MPI](https://www.open-mpi.org/) for parallelization
- [GNU Scientific Library](https://www.gnu.org/software/gsl/)
- [BLAS](https://www.openblas.net/)
- C++-17

### Visualization
- [gnuplot](http://gnuplot.info)
- [gnuplot-iostream](https://github.com/dstahlke/gnuplot-iostream)
- [Boost](https://www.boost.org/)
  - [Boost.Iostreams](https://www.boost.org/doc/libs/1_80_0/libs/iostreams/doc/index.html)
  - [Boost.Filesystem](https://www.boost.org/doc/libs/1_78_0/libs/filesystem/doc/index.htm)
- [Eigen](https://eigen.tuxfamily.org/)

## Configuration and Building
1. Clone the repository
2. Copy the file config.mk.in to config.mk
3. Edit the contents of config.mk to correspond to the locations on your system of the above dependencies

After navigating into the repository and creating/editing config.mk as described above, `make` will build the executables. If you're only interested in building the simulation code, you can run `make swim_steady` to build the solver without the visualization tools (avoiding the need to install their dependencies.) Similarly, `make gen_gp` will build only the visualization tools.

### Running
To solve a system with separation distance S, left sphere radius lr, right sphere radius rr, and inertial parameter M^2, run

`mpirun -np num_procs swim_steady -s S -R lr -r rr -M M^2 -o output_directory m p`

to solve the system over a 2m x m bispherical grid of finite elements. The velocity elements are of polynomial order p, and the pressure elements are of order p-1. The code will print to stdout:
- the provided M^2 value
- the density of the left sphere
- the density of the right sphere
- the resulting steady swim speed

Run `./swim_steady -h` for a short help message and `./swim_steady -H` for more documentation on option parameters. Within the output directory, the Brinkman and Stokes flow fields will also be saved as binary files. The executable `./gen_gp` can be used to read in these files and construct heat maps of the Stokes vorticity field with overlaid streamlines.

## Example
The command `mpirun -np 4 swim_steady -s 3 -R 0.5 -r 1 -M 10 -o out 64 3` produces the output

`10.00000000 10.00000000 10.00000000 -0.00483433`

and the command `./gen_gp -o figure -A 12 -I -P 800 -T 100 -d 0.001 -X 12 -Y 12 -q vel -s w -p r out_s.bin` produces the image ![figure.pdf](./figure.pdf), which is 800 pixels across and contains 12 columns of streamlines. Its horizontal and vertical coordinate widths are both 12. Color denotes vorticity of the steady Stokes flow field, which the streamlines are derived from. The streamlines are integrated for dimensionless time T=100, using a local error tolerance of 0.001 for adaptive timestepping (i.e. the lines become longer as `-T` takes larger values and smoother as `-d` takes smaller values.)

Note that on some machines it may be required to set `OMP_NUM_THREADS=1` when running `swim_steady`.

## Contact
For questions, contact [Nick Derr](mailto:derr@mit.edu?subject=[reciprocal_swimming]). 
