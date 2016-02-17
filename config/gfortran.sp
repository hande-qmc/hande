[DEFAULT]
include_f = -I $${HDF5_ROOT-/usr}/include

[main]
fc = gfortran
cc = gcc
cxx = g++
cppflags = -DSINGLE_PRECISION
ld = gfortran
ldflags = -L ${HOME}/lib -L $${HDF5_ROOT-/usr}/lib -L/usr/local/lib
libs = -ltrlan -llapack -lblas -lstdc++ -lhdf5_fortran -lhdf5 -lz -luuid -llua -ldl
f90_module_flag = -J

[opt]
fflags = %(include_f)s -O3 -mfpmath=sse -msse2
cflags = -O3 -mfpmath=sse -msse2
cxxflags = -O3 -mfpmath=sse -msse2

[dbg]
fflags = %(include_f)s -g -fbounds-check -Wall -Wextra -fbacktrace -mfpmath=sse -msse2
cflags = -g -Wall -Wextra -fbacktrace -mfpmath=sse -msse2
cxxflags = -g -Wall -Wextra -fbacktrace -mfpmath=sse -msse2
