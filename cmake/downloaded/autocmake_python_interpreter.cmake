# (c) https://github.com/coderefinery/autocmake/blob/master/AUTHORS.md
# licensed under BSD-3: https://github.com/coderefinery/autocmake/blob/master/LICENSE

#.rst:
#
# Detects Python interpreter.
#
# Variables used::
#
#   PYTHON_INTERPRETER         - User-set path to the Python interpreter
#
# Variables defined::
#
#   PYTHONINTERP_FOUND         - Was the Python executable found
#   PYTHON_EXECUTABLE          - path to the Python interpreter
#   PYTHON_VERSION_STRING      - Python version found e.g. 2.5.2
#   PYTHON_VERSION_MAJOR       - Python major version found e.g. 2
#   PYTHON_VERSION_MINOR       - Python minor version found e.g. 5
#   PYTHON_VERSION_PATCH       - Python patch version found e.g. 2
#
# autocmake.yml configuration::
#
#   docopt: "--python=<PYTHON_INTERPRETER> The Python interpreter (development version) to use. [default: '']."
#   define: "'-DPYTHON_INTERPRETER=\"{0}\"'.format(arguments['--python'])"

if("${PYTHON_INTERPRETER}" STREQUAL "")
    find_package(PythonInterp REQUIRED)
else()
    if(NOT EXISTS "${PYTHON_INTERPRETER}")
        find_program(PYTHON_EXECUTABLE NAMES ${PYTHON_INTERPRETER})
        if (NOT EXISTS "${PYTHON_EXECUTABLE}")
            set(PYTHONINTERP_FOUND FALSE)
        endif()
    else()
        set(PYTHONINTERP_FOUND TRUE)
        set(PYTHON_EXECUTABLE "${PYTHON_INTERPRETER}")
    endif()
endif()
find_package(PythonInterp REQUIRED)
