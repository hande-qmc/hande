#.rst:
#
# Creates print_info.c and git_info.f90 in the build directory.
#

function(generate_info_header _location _name)
  # _location: where the info file should be generated
  # _name: the info file name, complete with extension
  # Accepts on optional argument for the name of the Fortran module holding the Git version info
  find_package(Git QUIET)

  set(_git_last_commit_hash "unknown")
  set(_git_last_commit_author "unknown")
  set(_git_last_commit_date "unknown")
  set(_git_branch "unknown")
  set(_git_describe "1.1-dev")

  if(GIT_FOUND)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} --no-pager show -s --pretty=format:%H -n 1
      OUTPUT_VARIABLE _git_last_commit_hash
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
      )

    execute_process(
      COMMAND ${GIT_EXECUTABLE} --no-pager show -s --pretty=format:%aN -n 1
      OUTPUT_VARIABLE _git_last_commit_author
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
      )

    execute_process(
      COMMAND ${GIT_EXECUTABLE} --no-pager show -s --pretty=format:%ad -n 1
      OUTPUT_VARIABLE _git_last_commit_date
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
      )

    execute_process(
      COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
      OUTPUT_VARIABLE _git_branch
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
      )

    execute_process(
      COMMAND ${GIT_EXECUTABLE} describe --abbrev=7 --long --always --dirty --tags
      OUTPUT_VARIABLE _git_describe
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
      )
  endif()

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
    set(_pop_size ${HANDE_POP_SIZE})
  endif()

  set(_det_size "unknown")
  if(DET_SIZE)
    set(_det_size ${HANDE_DET_SIZE})
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
  set(_use_uuid "DISABLE_UUID defined.  UUID disabled.")
  if(USE_UUID)
    set(_use_uuid "DISABLE_UUID not defined.  UUID enabled.")
  endif()
  # ScaLAPACK
  set(_use_scalapack "DISABLE_SCALAPACK defined.  ScaLAPACK disabled.")
  if(USE_ScaLAPACK)
    set(_use_scalapack "DISABLE_SCALAPACK not defined.  ScaLAPACK enabled.")
  endif()
  # SINGLE_PRECISION
  set(_use_single_precision "SINGLE_PRECISION not defined.  Double precision used throughout.")
  if(USE_SINGLE_PRECISION)
    set(_use_single_precision "SINGLE_PRECISION defined.  Single precision used where relevant.")
  endif()
  # POP_CNT
  set(_use_popcnt "USE_POPCNT not defined.  Internal POPCNT procedure used.")
  if(USE_INTRINSIC_POPCNT)
    set(_use_popcnt "USE_POPCNT defined.  Fortran 2003 POPCNT procedure used.")
  endif()
  # dSFMT Mersenne Exponent
  set(_dsfmt_mexp ${HANDE_DSFMT_MEXP})

  # Parallel setup
  set(_enable_mpi ${USE_MPI})
  set(_mpi_launcher "unknown")
  if(USE_MPI)
    set(_mpi_launcher ${MPIEXEC})
  endif()
  set(_enable_omp ${USE_OPENMP})

  configure_file(
    ${PROJECT_SOURCE_DIR}/cmake/custom/print_info.c.in
    ${_location}/${_name}
    @ONLY
    )

  # Check if the optional argument was passed
  if(ARGN)
    configure_file(
      ${PROJECT_SOURCE_DIR}/cmake/custom/git_info.f90.in
      ${_location}/${ARGV2}
      @ONLY
      )
  endif()

  unset(_git_last_commit_hash)
  unset(_git_last_commit_author)
  unset(_git_last_commit_date)
  unset(_git_branch)
  unset(_git_describe)

  add_custom_target(
    build_info
    ALL DEPENDS ${_location}/${_name}
    )
endfunction()
