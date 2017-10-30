#include <stdio.h>
#include <unistd.h>

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
  printf("Compilation hostname     | %s\n", HOST_NAME);
  printf("Operating system         | %s\n", SYSTEM);
  printf("CMake version            | %s\n", CMAKE_VERSION);
  printf("CMake generator          | %s\n", CMAKE_GENERATOR);
  printf("CMake build type         | %s\n", CMAKE_BUILD_TYPE);
  printf("Configuration time       | %s\n", CONFIGURATION_TIME);
  printf("Python version           | %s\n", PYTHON_VERSION);
  printf("Fortran compiler         | %s\n", FORTRAN_COMPILER);
  printf("C compiler               | %s\n", C_COMPILER);
  printf("C++ compiler             | %s\n", CXX_COMPILER);
  printf("DET_SIZE set to          | %s\n", DET_SIZE);
  printf("POP_SIZE set to          | %s\n", POP_SIZE);
  printf("dSFMT Mersenne exponent  | %s\n", dSFMT_MEXP);
  printf("MPI parallelization      | %s\n", ENABLE_MPI);
  printf("MPI launcher             | %s\n", MPI_LAUNCHER);

  printf("\n");
  printf("\n");
  printf("Further components\n");
  printf("------------------\n");
  printf("\n");
  printf("   %s\n", HDF5);
  printf("   %s\n", LANCZOS);
  printf("   %s\n", UUID);
  printf("   %s\n", SCALAPACK);
  printf("   %s\n", POPCNT);
  printf("   %s\n", SINGLE_PRECISION);

  printf("\n");
  printf("\n");
  printf("Runtime information\n");
  printf("-------------------\n");
  printf("\n");
  char hostname[1024];
  int stat = gethostname(hostname, 1024);
  if (stat == 0) {
    printf("Hostname:                               \n    %s\n", hostname);
  } else {
    printf("Hostname (truncated to 1024 characters):\n    %s\n", hostname);
  }
  char cwd[1024];
  if (getcwd(cwd, sizeof(cwd)) != NULL) {
    printf("Current working dir:                    \n    %s\n", cwd);
  } else {
    perror("getcwd() error");
  }

  fflush(stdout);
}
