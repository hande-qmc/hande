#!/usr/bin/env bash

set -euo pipefail

cd "$HOME"/Downloads

echo "-- Installing TRLan"
if [[ ! -f $DEPS/trlan/lib/libtrlan.a ]]; then
  echo "-- TRLan NOT FOUND in cache"
  curl -Ls https://codeforge.lbl.gov/frs/download.php/210/trlan-201009.tar.gz | tar -xz
  (
   cd trlan-201009
   # Patch Make.inc to use gfortran-8
   patch < "$TRAVIS_BUILD_DIR"/.ci/Make.inc.patch
   make lib &> /dev/null
   # Install by hand
   mkdir -p "$DEPS"/trlan/lib
   cp libtrlan.a "$DEPS"/trlan/lib
  )
else
  echo "-- TRLan FOUND in cache"
fi
echo "-- Done with TRLan"

cd "$TRAVIS_BUILD_DIR"
