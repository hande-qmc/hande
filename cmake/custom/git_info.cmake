#.rst:
#
# Creates git_info.h in the build directory.
# This file can be included in Fortran sources to print
# Git repository version and status information
# to the program output.
#

function(generate_git_info_header _header_location _header_name)
  # _header_location: where the Git info header file should be generated
  # _header_name: the Git info header name, complete with extension (.h, .hpp, .hxx or whatever)
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

  configure_file(
    ${PROJECT_SOURCE_DIR}/cmake/custom/git_info.h.in
    ${_header_location}/${_header_name}
    @ONLY
    )

  # Check if the optional argument was passed
  if(ARGN)
    configure_file(
      ${PROJECT_SOURCE_DIR}/cmake/custom/git_info.f90.in
      ${_header_location}/${ARGV2}
      @ONLY
      )
  endif()

  unset(_git_last_commit_hash)
  unset(_git_last_commit_author)
  unset(_git_last_commit_date)
  unset(_git_branch)
  unset(_git_describe)

  add_custom_target(
    git_info
    ALL DEPENDS ${_header_location}/${_header_name}
    )
endfunction()
