def parse_cmake_module(s_in, overrides={}):
    import sys
    from collections import Mapping, Iterable, defaultdict
    from autocmake.parse_yaml import parse_yaml

    # we do not use the nicer sys.version_info.major
    # for compatibility with Python < 2.7
    if sys.version_info[0] > 2:
        from io import StringIO
    else:
        from StringIO import StringIO

    parsed_config = defaultdict(lambda: None)

    if 'autocmake.yml configuration::' not in s_in:
        return parsed_config

    s_out = []
    is_rst_line = False
    for line in s_in.split('\n'):
        if is_rst_line:
            if len(line) > 0:
                if line[0] != '#':
                    is_rst_line = False
            else:
                is_rst_line = False
        if is_rst_line:
            s_out.append(line[2:])
        if '#.rst:' in line:
            is_rst_line = True

    autocmake_entry = '\n'.join(s_out).split('autocmake.yml configuration::')[1]
    autocmake_entry = autocmake_entry.replace('\n  ', '\n')

    buf = StringIO(autocmake_entry)
    config = parse_yaml(buf, overrides)

    for k, v in config.items():
        if isinstance(v, Iterable) and not isinstance(v, str):
            parsed_config[k] = [x for x in config[k]]
        else:
            parsed_config[k] = [config[k]]

    return parsed_config


def test_parse_cmake_module():

    s = r'''#.rst:
#
# Foo ...
#
# autocmake.yml configuration::
#
#   docopt:
#     - "--cxx=<CXX> C++ compiler [default: g++]."
#     - "--extra-cxx-flags=<EXTRA_CXXFLAGS> Extra C++ compiler flags [default: '']."
#   export: "'CXX={0}'.format(arguments['--cxx'])"
#   define: "'-DEXTRA_CXXFLAGS=\"{0}\"'.format(arguments['--extra-cxx-flags'])"

enable_language(CXX)

if(NOT DEFINED CMAKE_C_COMPILER_ID)
    message(FATAL_ERROR "CMAKE_C_COMPILER_ID variable is not defined!")
endif()'''

    parsed_config = parse_cmake_module(s)
    assert parsed_config['docopt'] == ["--cxx=<CXX> C++ compiler [default: g++].", "--extra-cxx-flags=<EXTRA_CXXFLAGS> Extra C++ compiler flags [default: '']."]


def test_parse_cmake_module_no_key():

    s = '''#.rst:
#
# Foo ...
#
# Bar ...

enable_language(CXX)

if(NOT DEFINED CMAKE_C_COMPILER_ID)
    message(FATAL_ERROR "CMAKE_C_COMPILER_ID variable is not defined!")
endif()'''

    parsed_config = parse_cmake_module(s)
    assert parsed_config['docopt'] is None


def test_parse_cmake_module_interpolate():

    s = r'''#.rst:
#
# Foo ...
#
# autocmake.yml configuration::
#
#   major: 1
#   minor: 2
#   patch: 3
#   a: v%(major)
#   b: v%(minor)
#   c: v%(patch)

enable_language(CXX)'''

    parsed_config = parse_cmake_module(s)
    assert parsed_config['a'] == ['v1']
    assert parsed_config['b'] == ['v2']
    assert parsed_config['c'] == ['v3']


def test_parse_cmake_module_overrides():

    s = r'''#.rst:
#
# Foo ...
#
# autocmake.yml configuration::
#
#   major: 1
#   minor: 2
#   patch: 3
#   a: v%(major)
#   b: v%(minor)
#   c: v%(patch)

enable_language(CXX)'''

    parsed_config = parse_cmake_module(s, {'minor': 4})
    assert parsed_config['a'] == ['v1']
    assert parsed_config['b'] == ['v4']
    assert parsed_config['c'] == ['v3']
