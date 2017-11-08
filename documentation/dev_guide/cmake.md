# Rewriting HANDE's build system

This document is still a work-in-progress!!!

The contents of `src` are compiled into a library `libhande` and an executable `hande.x`.

I have dropped:
  1. `aotus-5.2` and Lua 5.2 support to make the transition to CMake easier.

## CMake options: user-facing and internals

Conventions:
- `ENABLE_<component>` are _boolean_, **user-facing** options to enable optional components.
  The CMake code will check for availability of additional dependencies
  and set a corresponding `USE_<component>` **internal** option to `ON` or `OFF`.
  These options can be mapped to preprocessor flags `-DENABLE_<component>`. These can be injected
  into the build system generator by using `target_compile_definitions` using the following
  pattern exploiting CMake generator expressions:
  ```
  target_compile_definitions(<target-name>
    PRIVATE
      $<$<BOOL:${USE_<option>}>:<option-preprocessor-variable>>
      $<$<NOT:$<BOOL:${USE_<option>>}>>:<option-preprocessor-variable>>
    )
  ```
  Options in this class:
    * `ENABLE_BACKTRACE`. Sets the `DISABLE_BACKTRACE` preprocessor variable.
    * `ENABLE_HDF5`. Sets the `DISABLE_HDF5` preprocessor variable.
    * `ENABLE_SINGLE_PRECISION`. Sets the `SINGLE_PRECISION` preprocessor variable.
    * `ENABLE_MPI`. Sets the `PARALLEL` preprocessor variable.
    * `ENABLE_OPENMP`.
    * `ENABLE_INTRINSIC_POPCNT`. Set the `USE_POPCNT` preprocessor variable.
    * `ENABLE_UUID`. Sets the `DISABLE_UUID` preprocessor variable.
    * `ENABLE_LANCZOS`. Sets the `DISABLE_LANCZOS` preprocessor variable.
    * `ENABLE_ScaLAPACK`. Sets the `DISABLE_SCALAPACK` preprocessor variable.
- `HANDE_<option>` are _non-boolean_, **user-facing** options. CMake will check for consistency
  against their available valid values. These options can be mapped to preprocessor `-D<option>` flags
  that can be injected into the build system generator by using `target_compile_definitons`:
  ```
  target_compile_definitions(<target-name>
    PRIVATE
      <option-preprocessor-variable>=${HANDE_<option>}
    )
  ```
  Options in this class:
    * `HANDE_DET_SIZE`. Its value is used to set the `DET_SIZE` preprocessor variable.
    * `HANDE_POP_SIZE`. Its value is used to set the `POP_SIZE` preprocessor varaible.
    * `HANDE_DSFMT_MEXP`. Its value is used to set the `DSFMT_MEXP` preprocessor variable.

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
