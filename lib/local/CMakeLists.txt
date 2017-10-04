list(APPEND sources
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

add_library(coreutils
  OBJECT
    ${sources}
  )
if(NOT USE_HDF5)
  target_compile_definitions(coreutils
    PUBLIC
      DISABLE_HDF5
    )
endif()
set_target_properties(coreutils
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
  )
add_library(HANDE::coreutils ALIAS coreutils)
