#include <stdio.h>

#include "build_info.h"
#include "git_info.h"

/*
 * Print the configuration/build info generated at build time by the executable.
 * We do this from C and not from Fortran because Fortran has trouble with long
 * lines.
 */

void print_info(void) {
  printf("\n");
  printf("\n");
  printf("Version information\n");
  printf("-------------------\n");
  printf("\n");
  printf("Version       | %s\n", GIT_DESCRIBE);
  printf("Commit hash   | %s\n", GIT_COMMIT_HASH);
  printf("Commit author | %s\n", GIT_COMMIT_AUTHOR);
  printf("Commit date   | %s\n", GIT_COMMIT_DATE);
  printf("Branch        | %s\n", GIT_BRANCH);

  printf("\n");
  printf("\n");
  printf("Configuration and build information\n");
  printf("-----------------------------------\n");
  printf("\n");
  printf("Who compiled             | %s\n", USER_NAME);
  printf("Compiled on server       | %s\n", HOST_NAME);
  printf("Operating system         | %s\n", SYSTEM);
  printf("CMake version            | %s\n", CMAKE_VERSION);
  printf("CMake generator          | %s\n", CMAKE_GENERATOR);
  printf("CMake build type         | %s\n", CMAKE_BUILD_TYPE);
  printf("Configuration time       | %s\n", CONFIGURATION_TIME);
  printf("Python version           | %s\n", PYTHON_VERSION);
  printf("Fortran compiler         | %s\n", FORTRAN_COMPILER);
  printf("C compiler               | %s\n", C_COMPILER);
  printf("C++ compiler             | %s\n", CXX_COMPILER);
  printf("MPI parallelization      | %s\n", ENABLE_MPI);
  printf("MPI launcher             | %s\n", MPI_LAUNCHER);
  printf("\n");
  printf("\n");

  fflush(stdout);
}
