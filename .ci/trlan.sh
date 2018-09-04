#!/usr/bin/env bash

set -euo pipefail

cd "$HOME"/Downloads

echo "-- Installing TRLan"
if [[ ! -f $DEPS/trlan/lib/libtrlan.a ]]; then
  echo "-- TRLan NOT FOUND in cache"
  curl -Ls https://codeforge.lbl.gov/frs/download.php/210/trlan-201009.tar.gz | tar -xz
  (
   cd trlan-201009
   env CC=gcc-8 FC=gfortran-8 make lib &> /dev/null
   # Install by hand
   mkdir -p "$DEPS"/trlan/lib
   cp libtrlan.a "$DEPS"/trlan/lib
  )
else
  echo "-- TRLan FOUND in cache"
fi
echo "-- Done with TRLan"

cd "$TRAVIS_BUILD_DIR"
