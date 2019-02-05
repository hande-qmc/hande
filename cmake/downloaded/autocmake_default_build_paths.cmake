# (c) https://github.com/dev-cafe/autocmake/blob/master/AUTHORS.md
# licensed under BSD-3: https://github.com/dev-cafe/autocmake/blob/master/LICENSE

#.rst:
#
# Sets binary and library output directories to ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
# and ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}, respectively.
#
# Variables modified::
#
#   CMAKE_ARCHIVE_OUTPUT_DIRECTORY
#   CMAKE_LIBRARY_OUTPUT_DIRECTORY
#   CMAKE_RUNTIME_OUTPUT_DIRECTORY

include(GNUInstallDirs)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
