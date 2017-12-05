# (c) https://github.com/coderefinery/autocmake/blob/master/AUTHORS.md
# licensed under BSD-3: https://github.com/coderefinery/autocmake/blob/master/LICENSE

#.rst:
#
# Sets binary and library output directories to ${PROJECT_BINARY_DIR}/bin
# and ${PROJECT_BINARY_DIR}/lib, respectively.
#
# Variables modified::
#
#   CMAKE_ARCHIVE_OUTPUT_DIRECTORY
#   CMAKE_LIBRARY_OUTPUT_DIRECTORY
#   CMAKE_RUNTIME_OUTPUT_DIRECTORY

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
