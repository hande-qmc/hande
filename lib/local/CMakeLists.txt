generate_info_header()

add_library(coreutils OBJECT
  alias.f90
  base_types.f90
  bit_table_256.f90
  checking.f90
  const.F90
  csr.f90
  environment_report.F90
  error_handling.F90
  hash.f90
  hash_table.f90
  hdf5_helper.F90
  json_out.f90
  linalg.F90
  parallel.F90
  ranking.f90
  shmem.F90
  utils.F90
  )
target_include_directories(coreutils
  PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}
  PUBLIC
    $<$<BOOL:${USE_HDF5}>:${HDF5_INCLUDE_DIRS}>
  )
target_sources(coreutils
  PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/git_info.f90
    ${CMAKE_CURRENT_BINARY_DIR}/print_info.c
  )
# We need to pass all compiler definitions to have a coherent set up report in
# print_info.c
target_compile_definitions(coreutils
  PRIVATE
    DET_SIZE=${HANDE_DET_SIZE}
    POP_SIZE=${HANDE_POP_SIZE}
    DSFMT_MEXP=${HANDE_DSFMT_MEXP}
    $<$<NOT:$<BOOL:${USE_BACKTRACE}>>:DISABLE_BACKTRACE>
    $<$<NOT:$<BOOL:${USE_HDF5}>>:DISABLE_HDF5>
    $<$<NOT:$<BOOL:${USE_UUID}>>:DISABLE_UUID>
    $<$<NOT:$<BOOL:${USE_ScaLAPACK}>>:DISABLE_SCALAPACK>
    $<$<BOOL:${USE_INTRINSIC_POPCNT}>:USE_POPCNT>
    $<$<BOOL:${USE_SINGLE_PRECISION}>:SINGLE_PRECISION>
    $<$<BOOL:${USE_MPI}>:PARALLEL>
    $<$<NOT:$<BOOL:${HAS_MPI3}>>:DISABLE_MPI3>
    $<$<CONFIG:Debug>:DEBUG>
  )
add_library(HANDE::coreutils ALIAS coreutils)
