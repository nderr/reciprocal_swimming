###############################################################################
# Active mixtures simulation makefile
#
# Rycroft Group, Harvard SEAS
#
# Author : Nick Derr
# Date   : Aug 1, 2019
#
# Note: The include statements below generate the cflags, iflags, and lflags
# variables used in the compilation call. You may have to monkey with this
# in the case of circular dependencies.
#
# Targets: all executables libraries depend clean
###############################################################################

#######################################
# knobs to turn!

# 1 on, 0 off
DEBUG:=0

# 0 off, 1-3 on
OPT:=3

# 0 complex, 1 real
REAL:=0
#######################################

# executables
solves:=swim_steady
stats:=gen_gp
execs:=$(solves) $(stats)

# libraries comprising this project
libs:=stats args stokes

# objects comprising each library
stats_objs=gp_data.o
args_objs=arg_helpers.o col_writer.o arg_manager.o
stokes_objs=fem.o quad.o poly.o master.o lapack.o brinkman.o

# flags required for each library
stats_flags=
args_flags=
stokes_flags=-lgsl

#####################################

# put together base flags
cflags=-std=c++17 -fopenmp
iflags:=
lflags:=-L.
xflags:=

#####################################

include config.mk

# are we debugging?
ifeq ($(strip $(DEBUG)),1)
	petsc_arch_complex:=$(petsc_arch_complex_debug)
else
	petsc_arch_complex:=$(petsc_arch_complex_opt)
endif

# define locations of real and complex petsc installs

# include this info in i,l,x flags
iflags += -I$(petsc_dir)/include
iflags += -I$(petsc_arch_complex)/include
lflags += -L$(petsc_arch_complex)/lib
xflags += -Wl,-rpath $(petsc_arch_complex)/lib

# add gnuplot
iflags += -I$(gnuplot_io_loc)

# add gsl
iflags += -I$(gsl_I)
lflags += -L$(gls_L)

# add blas
lflags += -L$(blas_L) $(blas_lflag)

# optimization/debugging values
ifeq ($(strip $(DEBUG)),1)
	dflag=-g
endif
oflag=-O$(OPT)

# pull in petsc locations, build libraries
include Makefile.petsc
include Makefile.in

# rules for each exectuable
flags=$(cflags) $(dflag) $(oflag) $(iflags) $(lflags) $(xflags)
compile=$(mpicxx) $(flags) -o $@ $< 
lfem=-lgsl -lgslcblas -llapack
base=-largs -lstokes $(lfem)

$(solves): %: %.cc libargs.a libstokes.a
	$(compile) $(base) -lpetsc -llapack

$(stats): %: %.cc libargs.a libstokes.a libstats.a
	$(compile) $(base) -lstats -lpetsc -lboost_iostreams -lboost_filesystem -llapack
