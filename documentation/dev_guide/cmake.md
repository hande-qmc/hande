# Rewriting HANDE's build system

## Directory structure

`lib` has been split into two subdirectories:
  - `lib`. Contains all the source files that were previously under `lib/local`
  - `external`. Contains all the sources that were previously under `lib` except stuff in `lib/local`:
    1. `aotus`.
    2. `dSFMT-src-2.2`
    3. `dSFMT_F03_interface`
    4. `external` This directory has been renamed to... _find a name_
    5. `cephes`
  I have dropped:
    1. `aotus-5.2` and Lua 5.2 support to make the transition to CMake easier.

The contents of `src` are compiled into a library `libhande` and an executable `hande.x`.
The directory is split into subdirectories, to give the code some more structure.

## Structure of src

I have tried to organize source files into a hierarchical structure that should
help avoid circular dependencies. All subdirectories but `core` may depend on stuff
in other subdirectories in `src`. The dependencies are explicitly listed in the `CMakeLists.txt`.
- `core` contains source code that depends on stuff in `lib`
- `top_level` contains source code that doesn’t clearly belong to any other
  subdirectory _and_ uses a lot of stuff from many other modules.
- `qmc`
- `dmqmc`
- `fciqmc`
- `ccmc`
- `top_level`

## Files renamed

- `core.f90` has become `hande.f90`
- `ccmc_data.F90` has become `multispawn.F90`. The `bit_string_ptr` and
  `cluster_t` types have been moved elsewhere.
- `qmc_io.f90` has become `hande_io.f90`

## Preprocessor definitions used in HANDE

- `POP_SIZE`
- `DET_SIZE`
- `PARALLEL`
- `SINGLE_PRECISION`
- `DEBUG`
- `DISABLE_BACKTRACE`
- `DISABLE_HDF5`
- `USE_POPCNT`
- `DISABLE_LANCZOS`

## Code changes

- Trivial: remove extra spaces.
- Trivial: fix typos in in-code comments.
- Types `bit_string_prt` and `cluster_t` have been moved to `determinants.f90`.
  This breaks the dependency between `determinants` and `ccmc_data`.
- Removed `dealloc.F90` by having the deallocation code where the data structure is defined.
