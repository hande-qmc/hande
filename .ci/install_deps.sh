#!/usr/bin/env bash

set -euo pipefail

echo "-- Download testcode"
git clone git://github.com/jsspencer/testcode.git $HOME/Deps/testcode &> /dev/null
echo "-- Done with testcode"

# OS-dependent operations
if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
  LUA_OS_NAME="linux"
elif [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
  LUA_OS_NAME="macosx"
  brew update &> /dev/null
  brew install gcc@7 ossp-uuid open-mpi cmake pyenv-virtualenv
fi

cd $HOME/Downloads

LUA_VERSION="5.3.4"
echo "-- Installing Lua $LUA_VERSION"
if [[ ! -f $LUA_ROOT/lib/liblua.a ]]; then
  echo "-- Lua $LUA_VERSION NOT FOUND in cache"
  curl -Ls http://www.lua.org/ftp/lua-$LUA_VERSION.tar.gz | tar -xz
  cd lua-$LUA_VERSION
  make $LUA_OS_NAME install INSTALL_TOP=$LUA_ROOT &> /dev/null
  cd $HOME/Downloads
else
  echo "-- Lua $LUA_VERSION FOUND in cache"
fi
echo "-- Done with Lua $LUA_VERSION"

HDF5_VERSION="1.10.1"
echo "-- Installing HDF5 $HDF5_VERSION"
HDF5_URL="https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${HDF5_VERSION%.*}/hdf5-$HDF5_VERSION/src/hdf5-$HDF5_VERSION.tar.bz2"
if [[ ! -f $HDF5_ROOT/lib/libhdf5_fortran.so ]]; then
  echo "-- HDF5 $HDF5_VERSION NOT FOUND in cache"
  curl -Ls $HDF5_URL | tar -xj
  cd hdf5-$HDF5_VERSION
  ./configure --prefix=$HDF5_ROOT --enable-fortran &> /dev/null
  make install --jobs=2 &> /dev/null
  cd $HOME/Downloads
else
  echo "-- HDF5 $HDF5_VERSION FOUND in cache"
fi
echo "-- Done with HDF5 $HDF5_VERSION"

ScaLAPACK_VERSION="2.0.2"
echo "-- Installing ScaLAPACK $ScaLAPACK_VERSION"
if [[ ! -f $ScaLAPACK_ROOT/lib/libscalapack.a ]]; then
  echo "-- ScaLAPACK $ScaLAPACK_VERSION NOT FOUND in cache"
  curl -Ls http://www.netlib.org/scalapack/scalapack-$ScaLAPACK_VERSION.tgz | tar -xz
  cd scalapack-$ScaLAPACK_VERSION
  cmake -H. -Bbuild -DCMAKE_INSTALL_PREFIX=$ScaLAPACK_ROOT &> /dev/null
  cmake --build build -- install --jobs=2 &> /dev/null
  cd $HOME/Downloads
else
  echo "-- ScaLAPACK $ScaLAPACK_VERSION FOUND in cache"
fi
echo "-- Done with ScaLAPACK $ScaLAPACK_VERSION"

echo "-- Installing TRLan"
if [[ ! -f $DEPS/trlan/lib/libtrlan.a ]]; then
  echo "-- TRLan NOT FOUND in cache"
  curl -Ls https://codeforge.lbl.gov/frs/download.php/210/trlan-201009.tar.gz | tar -xz
  cd trlan-201009
  make lib &> /dev/null
  # Install by hand
  mkdir -p $DEPS/trlan/lib
  cp libtrlan.a $DEPS/trlan/lib
  cd $HOME/Downloads
else
  echo "-- TRLan FOUND in cache"
fi
echo "-- Done with TRLan"

cd $TRAVIS_BUILD_DIR
