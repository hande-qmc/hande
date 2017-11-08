generate_info_header(${CMAKE_CURRENT_BINARY_DIR} print_info.c git_info.f90)

add_library(coreutils
  OBJECT
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
    utils.F90
  )

target_include_directories(coreutils
  PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}
  )
target_sources(coreutils
  PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/git_info.f90
    ${CMAKE_CURRENT_BINARY_DIR}/print_info.c
  )
target_compile_definitions(coreutils
  PUBLIC
    DET_SIZE=${HANDE_DET_SIZE}
    POP_SIZE=${HANDE_POP_SIZE}
    $<$<NOT:$<BOOL:${USE_BACKTRACE}>>:DISABLE_BACKTRACE>
    $<$<NOT:$<BOOL:${USE_UUID}>>:DISABLE_UUID>
    $<$<NOT:$<BOOL:${USE_HDF5}>>:DISABLE_HDF5>
    $<$<BOOL:${USE_SINGLE_PRECISION}>:SINGLE_PRECISION>
  )
set_target_properties(coreutils
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
  )

add_library(HANDE::coreutils ALIAS coreutils)
