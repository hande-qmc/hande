#!/usr/bin/env bash

set -euo pipefail

# OS-dependent operations
if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
    LUA_OS_NAME="linux"
elif [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    LUA_OS_NAME="macosx"
fi

cd "$HOME"/Downloads

LUA_VERSION="5.3.4"
echo "-- Installing Lua $LUA_VERSION"
if [[ ! -f $LUA_ROOT/lib/liblua.a ]]; then
    echo "-- Lua $LUA_VERSION NOT FOUND in cache"
    curl -Ls http://www.lua.org/ftp/lua-$LUA_VERSION.tar.gz | tar -xz
    (
     cd lua-$LUA_VERSION
     env CC=gcc-8 make $LUA_OS_NAME install INSTALL_TOP="$LUA_ROOT" &> /dev/null
    )
else
    echo "-- Lua $LUA_VERSION FOUND in cache"
fi
echo "-- Done with Lua $LUA_VERSION"

cd "$TRAVIS_BUILD_DIR"
