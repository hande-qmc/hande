#.rst:
#
# Creates print_info.c and git_info.f90 in the build directory.
#

function(generate_info_header)
  set(_mpi_launcher "unknown")
  if(USE_MPI)
    set(_mpi_launcher ${MPIEXEC})
  endif()

  file(COPY ${PROJECT_SOURCE_DIR}/tools/configurator.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  set(_config_script
"
import os, sys
sys.path.append('${CMAKE_CURRENT_BINARY_DIR}')
import configurator as cf
conf_dict = cf.prepare_configuration_dictionary(cmake_version='${CMAKE_VERSION}',
            cmake_generator='${CMAKE_GENERATOR}',
            Fortran_compiler='${CMAKE_Fortran_COMPILER}',
            C_compiler='${CMAKE_C_COMPILER}',
            build_type='${CMAKE_BUILD_TYPE}',
            mpi_launcher='${_mpi_launcher}',
            lua_version='\"${LUA_VERSION_STRING}\"',
            hdf5_version='\"${HDF5_VERSION}\"')
cf.configure_file(conf_dict, 'print_info.c', in_path=os.path.join('${PROJECT_SOURCE_DIR}', 'lib/local'), suffix='.in')
cf.configure_file(conf_dict, 'git_info.f90', in_path=os.path.join('${PROJECT_SOURCE_DIR}', 'lib/local'), suffix='.in')
")

  execute_process(
    COMMAND
    ${PYTHON_EXECUTABLE} "-c" ${_config_script}
    )

  add_custom_target(
    build_info
    ALL DEPENDS
      ${PROJECT_SOURCE_DIR}/tools/configurator.py
      ${PROJECT_SOURCE_DIR}/lib/local/print_info.c
      ${PROJECT_SOURCE_DIR}/lib/local/git_info.f90
    )
endfunction()
