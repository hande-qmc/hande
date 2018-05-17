# (c) https://github.com/coderefinery/autocmake/blob/master/AUTHORS.md
# licensed under BSD-3: https://github.com/coderefinery/autocmake/blob/master/LICENSE

#.rst:
#
# Add preprocessor definitions (example: --add-definitions="-DTHIS -DTHAT=137").
# These are passed all the way down to the compiler.
#
# Variables used::
#
#   PREPROCESSOR_DEFINITIONS
#
# autocmake.yml configuration::
#
#   docopt: "--add-definitions=<STRING> Add preprocesor definitions [default: '']."
#   define: "'-DPREPROCESSOR_DEFINITIONS=\"{0}\"'.format(arguments['--add-definitions'])"

if(NOT "${PREPROCESSOR_DEFINITIONS}" STREQUAL "")
    add_definitions(${PREPROCESSOR_DEFINITIONS})
endif()
