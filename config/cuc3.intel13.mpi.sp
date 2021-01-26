[DEFAULT]
cppflags_opt = -DHAVE_SSE2 -DPARALLEL -DSINGLE_PRECISION

[main]
fc = mpif90
cc = mpicc
cxx = mpiCC
ld = mpif90
libs = -mkl -lmkl_scalapack_lp64 -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_openmpi_lp64 -Wl,--end-group -lpthread -cxxlib -lhdf5_fortran -luuid -llua
f90_module_flag = -module

[opt]
cppflags = %(cppflags_opt)s
fflags = -O3
cxxflags = -O3

[dbg]
cppflags = %(cppflags_opt)s -DDEBUG
fflags = -g -traceback -CB
cxxflags = -g -traceback 
