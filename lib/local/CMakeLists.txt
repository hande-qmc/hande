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
  json_out.f90
  linalg.F90
  parallel.F90
  ranking.f90
  utils.F90
  )

if(USE_HDF5)
  list(APPEND sources
    hdf5_helper.F90
    )
endif()

add_library(hande-utils
  OBJECT
    ${sources}
  )
set_target_properties(hande-utils
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
  )
