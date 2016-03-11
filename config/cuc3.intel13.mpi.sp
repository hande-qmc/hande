[main]
fc = mpif90
cc = mpicc
cxx = mpiCC
ld = mpif90
cppflags = -DHAVE_SSE2 -DPARALLEL -DSINGLE_PRECISION
libs = -mkl -lmkl_scalapack_lp64 -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_openmpi_lp64 -Wl,--end-group -lpthread -cxxlib -lhdf5_fortran -luuid -llua -ltrlan_mpi
f90_module_flag = -module

[opt]
fflags = -O3
cxxflags = -O3

[dbg]
fflags = -g -traceback -CB
cxxflags = -g -traceback 
