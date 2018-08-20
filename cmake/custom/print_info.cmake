#.rst:
#
# Creates print_info.c and git_info.f90 in the build directory.
#

function(generate_info_header)
  set(_mpi_launcher "unknown")
  if(USE_MPI)
    set(_mpi_launcher ${MPIEXEC})
  endif()

  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/config_script.py
"import os, sys
sys.path.append('${PROJECT_SOURCE_DIR}/tools')
import configurator as cf
conf_dict = cf.prepare_configuration_dictionary(cmake_version='${CMAKE_VERSION}',
            cmake_generator='${CMAKE_GENERATOR}', 
            Fortran_compiler='${CMAKE_Fortran_COMPILER}',
            C_compiler='${CMAKE_C_COMPILER}',
            build_type='${CMAKE_BUILD_TYPE}',
            mpi_launcher='${_mpi_launcher}',
            lua_version='${LUA_VERSION_STRING}',
hdf5_version='${HDF5_VERSION}')
cf.configure_file(conf_dict, 'print_info.c', in_path='${CMAKE_CURRENT_SOURCE_DIR}', out_path='${CMAKE_CURRENT_BINARY_DIR}', suffix='.in')
cf.configure_file(conf_dict, 'git_info.f90', in_path='${CMAKE_CURRENT_SOURCE_DIR}', out_path='${CMAKE_CURRENT_BINARY_DIR}', suffix='.in')
")

  add_custom_command(
    OUTPUT
      print_info.c
      git_info.f90
    COMMAND
      ${PYTHON_EXECUTABLE} config_script.py
    WORKING_DIRECTORY
      ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS
      ${CMAKE_CURRENT_SOURCE_DIR}/print_info.c.in
      ${CMAKE_CURRENT_SOURCE_DIR}/git_info.f90.in
    )
endfunction()
