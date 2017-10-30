#.rst:
#
# Creates the build_info.h file in the build directory.
# This file can be included into sources to print out
# build information variables to the executable program output.

function(generate_build_info_header _header_location _header_name)
  # _header_location: where the build info header file should be generated
  # _header_name: the build info header name, complete with extension (.h, .hpp, .hxx or whatever)
  set(_user_name "unknown")
  execute_process(
    COMMAND
      whoami
    TIMEOUT
      1
    OUTPUT_VARIABLE
      _user_name
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

  cmake_host_system_information(RESULT _host_name QUERY HOSTNAME)
  if(NOT _host_name)
    set(_host_name "unknown")
  endif()

  set(_system "unknown")
  if(CMAKE_SYSTEM)
      set(_system ${CMAKE_SYSTEM})
  endif()

  set(_cmake_version "unknown")
  if(CMAKE_VERSION)
      set(_cmake_version ${CMAKE_VERSION})
  endif()

  set(_cmake_generator "unknown")
  if(CMAKE_GENERATOR)
      set(_cmake_generator ${CMAKE_GENERATOR})
  endif()

  set(_cmake_build_type "unknown")
  if(CMAKE_BUILD_TYPE)
      set(_cmake_build_type ${CMAKE_BUILD_TYPE})
  endif()

  string(TIMESTAMP _configuration_time "%Y-%m-%d %H:%M:%S [UTC]" UTC)
  if(NOT _configuration_time)
    set(_configuration_time "unknown")
  endif()

  set(_python_version "unknown")
  if(PYTHONINTERP_FOUND)
    set(_python_version ${PYTHON_VERSION_STRING})
  endif()

  foreach(_lang Fortran C CXX)
    set(_${_lang}_compiler "unknown")
    if(CMAKE_${_lang}_COMPILER)
      set(_${_lang}_compiler ${CMAKE_${_lang}_COMPILER})
    endif()
  endforeach()

  set(_pop_size "unknown")
  if(POP_SIZE)
    set(_pop_size ${POP_SIZE})
  endif()

  set(_det_size "unknown")
  if(DET_SIZE)
    set(_det_size ${DET_SIZE})
  endif()

  # HDF5
  set(_use_hdf5 "DISABLE_HDF5 defined.  HDF5 disabled.")
  if(USE_HDF5)
    set(_use_hdf5 "DISABLE_HDF5 not defined.  HDF5 enabled.")
  endif()
  # LANCZOS
  set(_use_lanczos "DISABLE_LANCZOS defined.  Lanczos disabled.")
  if(USE_LANCZOS)
    set(_use_lanczos "DISABLE_LANCZOS not defined.  Lanczos enabled.")
  endif()
  # UUID
  set(_use_libuuid "DISABLE_UUID defined.  UUID disabled.")
  if(USE_LIBUUID)
    set(_use_libuuid "DISABLE_UUID not defined.  UUID enabled.")
  endif()
  # ScaLAPACK
  set(_use_scalapack "DISABLE_SCALAPACK defined.  ScaLAPACK disabled.")
  if(USE_SCALAPACK)
    set(_use_scalapack "DISABLE_SCALAPACK not defined.  ScaLAPACK enabled.")
  endif()
  # SINGLE_PRECISION
  set(_use_single_precision "SINGLE_PRECISION defined.  Single precision used where relevant.")
  if(USE_SINGLE_PRECISION)
    set(_use_single_precision "SINGLE_PRECISION not defined.  Double precision used throughout.")
  endif()
  # POP_CNT
  set(_use_popcnt "USE_POPCNT defined.  Fortran 2003 POPCNT procedure used.")
  if(USE_POPCNT)
    set(_use_popcnt "USE_POPCNT not defined.  Internal POPCNT procedure used.")
  endif()
  # dSFMT Mersenne Exponent
  set(_dsfmt_mexp ${DSFMT_MEXP})

  set(_enable_mpi ${ENABLE_MPI})

  set(_mpi_launcher "unknown")
  if(MPI_FOUND)
    set(_mpi_launcher ${MPIEXEC})
  endif()

  configure_file(
    ${PROJECT_SOURCE_DIR}/cmake/custom/build_info.h.in
    ${_header_location}/${_header_name}
    @ONLY
    )

  add_custom_target(
    build_info
    ALL DEPENDS ${_header_location}/${_header_name}
    )
endfunction()
