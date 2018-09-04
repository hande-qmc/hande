#!/usr/bin/env bash

set -euo pipefail

cd "$HOME"/Downloads

OpenMPI_VERSION="3.1.0"
echo "-- Installing OpenMPI $OpenMPI_VERSION"
OpenMPI_URL="https://download.open-mpi.org/release/open-mpi/v${OpenMPI_VERSION%.*}/openmpi-${OpenMPI_VERSION}.tar.gz"
if [[ ! -f $OpenMPI_ROOT/lib/libmpi.la ]]; then
    echo "-- OpenMPI $OpenMPI_VERSION NOT FOUND in cache"
    curl -Ls $OpenMPI_URL | tar -xz
    (
      cd openmpi-$OpenMPI_VERSION
      mkdir -p build
      cd build
      env CC=gcc-8 FC=gfortran-8 CXX=g++-8 ../configure --prefix="$OpenMPI_ROOT" &> /dev/null
      make all --jobs=2 &> /dev/null
      make install &> /dev/null
    )
else
    echo "-- OpenMPI $OpenMPI_VERSION FOUND in cache"
fi
echo "-- Done with OpenMPI $OpenMPI_VERSION"

cd "$TRAVIS_BUILD_DIR"
