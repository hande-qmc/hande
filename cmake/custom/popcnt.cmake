#.rst:
#
# Probes system for the intrinsic popcnt instruction.
#
# Variables used::
#
#   ENABLE_INTRINSIC_POPCNT
#
# Variables defined::
#
#   USE_INTRINSIC_POPCNT
#
# Variables modified::
#
#   CMAKE_Fortran_FLAGS
#
# autocmake.yml configuration::
#
#   docopt:
#     - "--popcnt Enable use of intrinsic popcnt [default: False]."
#   define: "'-DENABLE_INTRINSIC_POPCNT=\"{0}\"'.format(arguments['--popcnt'])"

option_with_print(ENABLE_INTRINSIC_POPCNT "Enable usage of popcnt intrinsic (requires hardware support)" OFF)
set(USE_INTRINSIC_POPCNT OFF)
if(ENABLE_INTRINSIC_POPCNT)
  set(_vendor_id)
  set(_cpu_family)
  set(_cpu_model)
  set(_cpu_flags)
  if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    file(READ "/proc/cpuinfo" _cpuinfo)
    string(REGEX REPLACE ".*vendor_id[ \t]*:[ \t]+([a-zA-Z0-9_-]+).*" "\\1" _vendor_id "${_cpuinfo}")
    string(REGEX REPLACE ".*cpu family[ \t]*:[ \t]+([a-zA-Z0-9_-]+).*" "\\1" _cpu_family "${_cpuinfo}")
    string(REGEX REPLACE ".*model[ \t]*:[ \t]+([a-zA-Z0-9_-]+).*" "\\1" _cpu_model "${_cpuinfo}")
    string(REGEX REPLACE ".*flags[ \t]*:[ \t]+([^\n]+).*" "\\1" _cpu_flags "${_cpuinfo}")
  elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    execute_process(COMMAND "/usr/sbin/sysctl -n machdep.cpu.vendor" OUTPUT_VARIABLE _vendor_id)
    execute_process(COMMAND "/usr/sbin/sysctl -n machdep.cpu.model"  OUTPUT_VARIABLE _cpu_model)
    execute_process(COMMAND "/usr/sbin/sysctl -n machdep.cpu.family" OUTPUT_VARIABLE _cpu_family)
    execute_process(COMMAND "/usr/sbin/sysctl -n machdep.cpu.features" OUTPUT_VARIABLE _cpu_flags)
    string(TOLOWER "${_cpu_flags}" _cpu_flags)
    string(REPLACE "." "_" _cpu_flags "${_cpu_flags}")
  endif()
  include(CheckFortranCompilerFlag)
  if(_cpu_flags MATCHES "popcnt")
    check_fortran_compiler_flag("-march=native" _march_native_fortran_works)
    check_fortran_compiler_flag("-xHost" _xhost_fortran_works)
    if(_march_native_fortran_works)
      set(USE_INTRINSIC_POPCNT ON)
      message(STATUS "Using processor's intrinsic popcnt instruction (-march=native compiler flag set)")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -march=native")
    elseif(_xhost_fortran_works)
      set(USE_INTRINSIC_POPCNT ON)
      message(STATUS "Using processor's intrinsic popcnt instruction (-xHost compiler flag set)")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -xHost")
    else()
      set(USE_INTRINSIC_POPCNT OFF)
      message(STATUS "No suitable Fortran compiler flag found: falling back to HANDE's own popcnt")
    endif()
  else()
    message(STATUS "Your processor doesn't supply popcnt: falling back to HANDE's own popcnt")
  endif()
endif()
