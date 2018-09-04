#!/usr/bin/env bash

set -euo pipefail

cd "$HOME"/Downloads

LAPACK_VERSION="3.8.0"
echo "-- Installing LAPACK $LAPACK_VERSION"
if [[ ! -f $LAPACK_ROOT/lib/liblapack.a ]]; then
  echo "-- LAPACK $LAPACK_VERSION NOT FOUND in cache"
  curl -Ls http://www.netlib.org/lapack/lapack-$LAPACK_VERSION.tar.gz | tar -xz
  (
   cd lapack-$LAPACK_VERSION
   cmake -H. -Bbuild -DCMAKE_INSTALL_PREFIX="$LAPACK_ROOT" -DCMAKE_C_COMPILER=gcc-8 -DCMAKE_Fortran_COMPILER=gfortran-8 &> /dev/null
   cmake --build build --target install -- --jobs=2 &> /dev/null
  )
else
  echo "-- LAPACK $LAPACK_VERSION FOUND in cache"
fi
echo "-- Done with LAPACK $LAPACK_VERSION"

cd "$TRAVIS_BUILD_DIR"
