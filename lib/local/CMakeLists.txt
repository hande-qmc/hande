generate_git_info_header(${CMAKE_CURRENT_BINARY_DIR} git_info.h git_info.f90)
generate_build_info_header(${CMAKE_CURRENT_BINARY_DIR} build_info.h)

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
    print_info.c
  )

target_include_directories(coreutils
  PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}
  )
target_sources(coreutils
  PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/git_info.f90
  )

set_target_properties(coreutils
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
  )

if(NOT USE_LIBUUID)
  target_compile_definitions(coreutils
    PUBLIC
      DISABLE_UUID
    )
endif()

if(NOT USE_HDF5)
  target_compile_definitions(coreutils
    PUBLIC
      DISABLE_HDF5
    )
endif()

add_library(HANDE::coreutils ALIAS coreutils)
