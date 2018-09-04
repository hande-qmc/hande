#!/usr/bin/env bash

set -euo pipefail

cd "$HOME"/Downloads

ScaLAPACK_VERSION="2.0.2"
echo "-- Installing ScaLAPACK $ScaLAPACK_VERSION"
if [[ ! -f $ScaLAPACK_ROOT/lib/libscalapack.a ]]; then
  echo "-- ScaLAPACK $ScaLAPACK_VERSION NOT FOUND in cache"
  curl -Ls http://www.netlib.org/scalapack/scalapack-$ScaLAPACK_VERSION.tgz | tar -xz
  (
   cd scalapack-$ScaLAPACK_VERSION
   cmake -H. -Bbuild -DCMAKE_INSTALL_PREFIX="$ScaLAPACK_ROOT" -DCMAKE_C_COMPILER=gcc-8 -DCMAKE_Fortran_COMPILER=gfortran-8 -DCMAKE_CXX_COMPILER=g++-8 -DMPI_BASE_DIR="$HOME"/Deps/openmpi &> /dev/null
   cmake --build build --target install -- --jobs=2 &> /dev/null
  )
else
  echo "-- ScaLAPACK $ScaLAPACK_VERSION FOUND in cache"
fi
echo "-- Done with ScaLAPACK $ScaLAPACK_VERSION"

cd "$TRAVIS_BUILD_DIR"
