#.rst:
#
# Sets exponent for the period of the Mersenne Twister (MT) random number
# generator (RNG)
#
# Variables modified::
#
#   DSFMT_MEXP
#
# autocmake.yml configuration::
#
#   docopt:
#     - "--dsfmt-mexp=<DSFMT_MEXP> An integer among 521, 1279, 2203, 4253, 11213,
#                        19937, 44497, 86243, 1322049, 216091 [default: 19937]."
#   define: "'-DDSFMT_MEXP=\"{0}\"'.format(arguments['--dsfmt-mexp'])"

message(STATUS "Exponent of the period of the Mersenne Twister RNG: ${DSFMT_MEXP}")
