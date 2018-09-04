#!/usr/bin/env bash

set -euo pipefail

cd "$HOME"/Downloads

HDF5_VERSION="1.10.1"
echo "-- Installing HDF5 $HDF5_VERSION"
HDF5_URL="https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${HDF5_VERSION%.*}/hdf5-$HDF5_VERSION/src/hdf5-$HDF5_VERSION.tar.bz2"
if [[ ! -f $HDF5_ROOT/lib/libhdf5_fortran.so ]]; then
    echo "-- HDF5 $HDF5_VERSION NOT FOUND in cache"
    curl -Ls $HDF5_URL | tar -xj
    (
     cd hdf5-$HDF5_VERSION
     env CC=gcc-8 FC=gfortran-8 CXX=g++-8 ./configure --prefix="$HDF5_ROOT" --enable-fortran &> /dev/null
     make install --jobs=2 &> /dev/null
    )
else
    echo "-- HDF5 $HDF5_VERSION FOUND in cache"
fi
echo "-- Done with HDF5 $HDF5_VERSION"

cd "$TRAVIS_BUILD_DIR"
