#!/usr/bin/env python
#
# Script for producing a make.inc file from a ini-style configuration file.
#
# http://github.com/jsspencer/generic_makefile
#
# copyright (c) 2012, James Spencer.
#
# MIT license:
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#
'''Produce a make.inc file for a specified target/configuration.

Usage: mkconfig [options] [config_file]

config_file can be either a path to a file or the name of a file in the
specified configuration directory.

If no configuration file is given then the script chooses a config file:
1) If ./make.inc exists and contains a line beginning
## mkconfig
then the remaining parameters on that line are used as argument.
Otherwise, if a .default file in the specified configuration directory
exists, it is used.

A configuration file does not need to be specified with the --ls option.

Options:
  -h, --help            show help message and exit.
  --help-long           show this help message and exit.
  -d CONFIG_DIR, --dir=CONFIG_DIR
                        Set directory containing the configuration files.
                        Default: config/.
  -g, --debug           Use the debug settings.  Default: use optimised
                        settings.
  -l, --ls              Print list of available configurations.
  -o OUT, --out=OUT     Set the output filename to which the makefile settings
                        are written.  Use -o - to write to stdout.  Default:
                        make.inc.

A platform is defined using a simple configuration file which is an ini file
consisting of three sections: main, opt and dbg.  For instance:

    [main]
    fc = gfortran
    ld = gfortran
    libs = -llapack -lblas

    [opt]
    fflags = -O3

    [dbg]
    fflags = -g

Any option not specified in the 'opt' and 'dbg' sections is inherited from the
'main' section.  The settings in 'opt' are used by default; the debug options
can be selected by passing the -g option to mkconfig.

All options are strings unless otherwise specified.  Available options are:

fc
    Set the fortran compiler.
fflags
    Set flags to be passed to the fortran compiler during compilation.
f90_module_flag
    Set the flag used by the compiler which is used to specify the directory
    where module (.mod) files are placed when created and where they should be
    searched for.
f90_module_flag_pad [boolean]
    True if a space needs to be inserted between the defined f90_module_flag
    and the corresponding directory argument.  Default: true.
cc
    Set the C compiler.
cflags
    Set flags to be passed to the C compiler during compilation.
ccd
    Set the C compiler used to generate the C dependency files.  Only required
    if cc doesn't support -MM and -MT flags.  Default: use cc.
cdflags
    Set the flags for the c++ compiler used to generate the C++ dependency files.
    Default: -MM -MT $CFLAGS
cxx
    Set the C++ compiler.
cxxflags
    Set flags to be passed to the C++ compiler during compilation.
cxxd
    Set the C compiler used to generate the C++ dependency files.  Only required
    if cc doesn't support -MM and -MT flags.  Default: use cxx.
cxxdflags
    Set the flags for the c++ compiler used to generate the C++ dependency files.
    Default: -MM -MT $CXXFLAGS
cpp
    Set the C preprocessor to be used on Fortran source files.  If not defined
    then the Fortran compiler is used to do the preprocessing.
cppflags
    Set flags to be used in the C preprocessing step.
    C preprocessing is applied to .F90, .F, .c and .cpp files (and not .f90
    files).
ld
    Set the linker program.
ldflags
    Set flags to be passed to the linker during linking of the compiled objects.
libs
    Set libraries to be used during the linking step.
ar
    Set the archive program.  Default: ar.
arflags
    Set the flags to be passed to the archive program.  Default: -rcs.
'''

try:
    import ConfigParser as configparser
except ImportError:
    import configparser
import optparse
import os
import sys
import shlex
import pipes
import contextlib
from configurator import prepare_configuration_dictionary, configure_file

MAKEFILE_TEMPLATE = '''# Generated by mkconfig.
## mkconfig %(args)s
# Using the %(config)s %(opt_level)s configuration.
#
# Not all options need to be set.

#-----
# Config info.

# Used to create a configuration-level and optimisation-level specific
# directory in which compiled objects are placed and to uniquely name binaries.

CONFIG = %(config)s
OPT = %(opt_level)s

LUA52 = %(lua_52)s

#-----
# Compiler configuration.

# --- Preprocessor ---
# If CPP is defined, then this is used to do the preprocessing on Fortran
# source files and then the Fortran compiler is called.  If CPP is not defined,
# then CPPFLAGS is passed to the Fortran compiler and the Fortran source files
# are preprocessed and compiled in a single step.
# All C and C++ files are preprocessed and compiled in one step by the C and
# C++ compilers, with CPPFLAGS passed to the C and C++ compilers.
CPP = %(cpp)s

# Preprocessing flags
# Used for *.F, *.F90, *.c and *.cpp files.
# Note three additional defintions specific to HANDE.
CPPFLAGS += %(cppflags)s

# --- Fortran ---
# compiler
FC = %(fc)s
# compiler flags
FFLAGS = %(fflags)s
# flag to specify directory to put compiled .mod files (e.g. -module for ifort, -J for gfortran)
F90_MOD_FLAG = %(f90_module_flag)s

# --- C ---
# compiler
CC = %(cc)s
# compiler flags
CFLAGS = %(cflags)s
# Not all compilers (e.g. xlc) support -MM -MT that we require to produce the
# C dependencies.  Fortunately there is (usually) a version of (e.g.) gcc
# kicking around which does.  If CCD is defined, then use that to generate the
# C dependency files, otherwise CC is used.
CCD = %(ccd)s
# If very custom flags are needed for the dependencies, they should be placed
# here, and will replace the -MM -MT $CFLAGS entirely
CCDFLAGS = %(ccdflags)s

# --- C++ ---
# compiler
CXX = %(cxx)s
# compiler flags
CXXFLAGS = %(cxxflags)s
# Not all compilers (e.g. xlc++) support -MM -MT that we require to produce the
# C++ dependencies.  Fortunately there is (usually) a version of (e.g.) g++
# kicking around which does.  If CXXD is defined, then use that to generate the
# C++ dependency files, otherwise CXX is used.
CXXD = %(cxxd)s
# If very custom flags are needed for the dependencies, they should be placed
# here, and will replace the -MM -MT $CFLAGS entirely
CXXDFLAGS = %(cxxdflags)s


# --- Linker ---
# linker
LD = %(ld)s
# linker flags
LDFLAGS = %(ldflags)s
# additional libraries to link to
LIBS = %(libs)s

# --- Archiver ---
# archive program (usually ar)
AR = %(ar)s
# archive flags (usually -rcs)
ARFLAGS = %(arflags)s
'''

def parse_options(my_args):
    '''Parse command line options given in the my_args list.

Returns the options object and the path to the selected configuration file.'''
    parser = optparse.OptionParser(usage='''mkconfig [options] [config_file]

config_file can be either a path to a file or the name of a file in the
specified configuration directory.

If no configuration file is given then the .default file in the specified
configuration directory is used.

A configuration file does not need to be specified with the --ls option.''')
    parser.add_option('--help-long', dest='help_long', action='store_true',
            help='Show long error message and exit.')
    parser.add_option('-d', '--dir', default='config/', dest='config_dir',
            help='Set directory containing the configuration files. '
                 'Default: %default.')
    parser.add_option('-g', '--debug', action='store_true', default=False,
            help='Use the debug settings.  Default: use optimised settings.')
    parser.add_option('-l', '--ls', action='store_true', default=False,
            help='Print list of available configurations.')
    parser.add_option('-o', '--out', default='make.inc',
            help='Set the output filename to which the makefile settings are '
                 'written.  Use -o - to write to stdout.  Default: %default.')
    parser.add_option('-D', '--define', action='append',
            help='Add e C preprocessor #define to be passed with -D to the compiler.'
                 ' Can be used multiple times.')

    (options, args) = parser.parse_args(my_args)

    # If there are no command line arguments, try and get them from make.inc
    if len(args) == 0 and os.path.exists('./make.inc'):
        with open('./make.inc') as f:
            for l in f:
                if '## mkconfig' in l:
                    my_args += shlex.split(l)[2:]
                    my_args = list(set(my_args))
                    print("Using arguments from:")
                    print(l.strip())
                    (options, args) = parser.parse_args(my_args)
                    break

    config_file = ''

    if not options.ls and not options.help_long:
        if len(args) == 1:
            config_file = args[0]
            if not os.path.exists(config_file):
                config_file = os.path.join(options.config_dir, config_file)
        elif len(args) == 0:
            if os.path.exists(os.path.join(options.config_dir, '.default')):
                config_file = os.path.join(options.config_dir, '.default')
            else:
                print('No config file provided.')
                parser.print_help()
                sys.exit(1)
        else:
            print('Too many config files provided.')
            parser.print_help()
            sys.exit(1)

    return (options, config_file, my_args)

def list_configs(config_dir):
    '''Return the path to all config files in config_dir.

We assume only config files are in config_dir.'''
    if os.path.isdir(config_dir):
        return [os.path.join(config_dir, f) for f in os.listdir(config_dir)]
    else:
        err = 'Config directory is not a directory: %s.' % (config_dir)
        raise IOError(err)

def parse_config(config_file):
    '''Parse the configuration file config_file.'''
    parser = configparser.SafeConfigParser()

    valid_sections = ['main', 'opt', 'dbg']

    valid_sections_upper = [section.upper() for section in valid_sections]

    valid_options = ['fc', 'fflags', 'cc', 'cflags', 'cxx', 'cxxflags', 'ccd',
            'cxxd', 'ccdflags', 'cxxdflags', 'cpp', 'cppflags', 'ld', 'ldflags', 'libs', 'ar', 'arflags',
            'f90_module_flag', 'f90_module_flag_pad', 'lua_52']

    boolean_states = {'0': False,
            '1': True,
            'false': False,
            'no': False,
            'off': False,
            'on': True,
            'true': True,
            'yes': True}

    default_settings = dict((key, '') for key in valid_options)
    default_settings.update(
            dict(ar='ar',
                 arflags='-rcs',
                 f90_module_flag_pad='True'
                )
            )

    parser.read(config_file)

    for section in parser.sections():
        if (section not in valid_sections and
                section not in valid_sections_upper):
            err = ('Invalid section in config file %s: %s.'
                    % (config_file, section))
            raise IOError(err)

    for section in valid_sections:
        if (section not in parser.sections() and
                section.upper() not in parser.sections()):
            err = ('Section not present in config file %s: %s.'
                        % (config_file, section))
            raise IOError(err)

    if not parser.sections():
        err = 'No sections in configuration file %s.' % (config_file)
        raise IOError(err)

    config = {}

    for section in ['opt', 'dbg']:
        config[section] = default_settings.copy()
        config[section].update(parser.items('main'))
        if section in parser.sections():
            config[section].update(parser.items(section))
        elif section.upper() in parser.sections():
            config[section].update(parser.items(section.upper()))
        # Update the module flag if a space needs to be appended.
        pad = config[section].pop('f90_module_flag_pad')
        # pad is a string rather than bool; use the same comparison that
        # configparser uses to determine if pad is set to true.
        if pad.lower() not in boolean_states:
            raise ValueError('Not a boolean: %s' % pad)
        if boolean_states[pad.lower()]:
            config[section]['f90_module_flag'] += ' '

    return config

def create_makefile(config_file, args, defs, use_debug=False):
    '''Returns the makefile using the options given in the config_file.

args is the list of arguments which is re-encoded in the make.inc.
    '''
    # Select the configuration section to use.
    if use_debug:
        config = parse_config(config_file)['dbg']
        config.update(opt_level='debug')
    else:
        config = parse_config(config_file)['opt']
        config.update(opt_level='optimised')

    # Set the name in the make.inc to be the basename of the config filename.
    # Follow any links.
    if os.path.islink(config_file):
        config_name = os.path.basename(os.readlink(config_file))
    else:
        config_name = os.path.basename(config_file)

    config.update(config=config_name)

    config["args"] = " ".join(pipes.quote(s) for s in args)
    if defs:
        config["defs"] = " -D".join([""]+defs)
    else:
        config["defs"] = ""
    conf_dict = prepare_configuration_dictionary(Fortran_compiler=config['fc'], build_type=config['opt_level'], C_compiler=config['cc'])
    configure_file(conf_dict, 'print_info.c', in_path=os.path.join(os.getcwd(), 'lib/local'), suffix='.in')
    configure_file(conf_dict, 'git_info.f90', in_path=os.path.join(os.getcwd(), 'lib/local'), suffix='.in')
    return (MAKEFILE_TEMPLATE % config)

@contextlib.contextmanager
def smart_open(filename=None):
    if filename and filename is not '-':
        fh = open(filename, 'w')
    else:
        fh = sys.stdout

    try:
        yield fh
    finally:
        if fh is not sys.stdout:
            fh.close()

def main(args):
    '''Main procedure.

args: list of arguments.'''

    (options, config_file, args) = parse_options(args)

    if options.help_long:
        print(__doc__)
    elif options.ls:
        configs = ', '.join(list_configs(options.config_dir))
        print('Available configurations in %s are: %s.'
                % (options.config_dir, configs))
    else:
        with smart_open(options.out) as fh:
            fh.write(create_makefile(config_file, args, options.define, options.debug))

if __name__ == '__main__':

    main(sys.argv[1:])
