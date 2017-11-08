#.rst:
#
# Sets exponent for the period of the Mersenne Twister (MT) random number
# generator (RNG)
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
