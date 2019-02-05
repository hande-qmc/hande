#.rst:
#
# Configure dSFMT library for use within HANDE::
#
#   - We set the exponent for the period of the Mersenne Twister (MT) random
#     number generator (RNG)
#   - We check whether the SSE2 instruction set is available and set the
#     HAVE_SSE2 flag accordingly.
#
# Variables modified::
#
#   HANDE_DSFMT_MEXP
#
# autocmake.yml configuration::
#
#   docopt:
#     - "--dsfmt-mexp=<HANDE_DSFMT_MEXP> An integer among 521, 1279, 2203, 4253, 11213,
#                        19937, 44497, 86243, 1322049, 216091 [default: 19937]."
#   define: "'-DHANDE_DSFMT_MEXP=\"{0}\"'.format(arguments['--dsfmt-mexp'])"

option_with_default(HANDE_DSFMT_MEXP "Exponent of the period of the Mersenne Twister RNG" 19937)
set(_VALID_DSFMT_MEXP 521 1279 2203 4253 11213 19937 44497 86243 1322049 216091)
if(DEFINED HANDE_DSFMT_MEXP AND NOT HANDE_DSFMT_MEXP IN_LIST _VALID_DSFMT_MEXP)
  message(STATUS "${HANDE_DSFMT_MEXP} not a valid exponent for a Mersenne prime, resetting to default value 19937")
  set(HANDE_DSFMT_MEXP 19937)
endif()

# Now check for the SSE2 instruction set by probing the system
function(system_has_sse2 _result)
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
     execute_process(
       COMMAND
         "/usr/sbin/sysctl" "-n" "machdep.cpu.vendor" "machdep.cpu.model" "machdep.cpu.family" "machdep.cpu.features"
       OUTPUT_VARIABLE
         _sysctl_output_string
       OUTPUT_STRIP_TRAILING_WHITESPACE
       )
     string(REPLACE "\n" ";" _sysctl_output ${_sysctl_output_string})
     list(GET _sysctl_output 0 _vendor_id)
     list(GET _sysctl_output 1 _cpu_model)
     list(GET _sysctl_output 2 _cpu_family)
     list(GET _sysctl_output 3 _cpu_flags)

     string(TOLOWER "${_cpu_flags}" _cpu_flags)
     string(REPLACE "." "_" _cpu_flags "${_cpu_flags}")
  endif()
  string(FIND _cpu_flags "sse2" _sse2_found)
  if(_sse2_found)
    set(${_result} TRUE PARENT_SCOPE)
    message(STATUS "CPU ${_vendor_id} with SSE2 instruction set FOUND")
  endif()
endfunction()