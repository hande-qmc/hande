#!/usr/bin/env bash

set -euo pipefail

cd "$HOME"/Downloads

CMAKE_VERSION="3.6.3"
echo "-- Installing CMake"
if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
  if [[ -f $HOME/Deps/cmake/$CMAKE_VERSION/bin/cmake ]]; then
    echo "-- CMake $CMAKE_VERSION FOUND in cache"
  else
    echo "-- CMake $CMAKE_VERSION NOT FOUND in cache"
    target_path=$HOME/Deps/cmake/$CMAKE_VERSION
    cmake_url="https://cmake.org/files/v${CMAKE_VERSION%.*}/cmake-${CMAKE_VERSION}-Linux-x86_64.tar.gz"
    mkdir -p "$target_path"
    curl -Ls $cmake_url | tar -xz -C "$target_path" --strip-components=1
  fi
elif [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
  brew upgrade cmake
fi
echo "-- Done installing CMake"

cd "$TRAVIS_BUILD_DIR"
