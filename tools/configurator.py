#!/usr/bin/env python
#
# Configures print_info.c.in and git_info.f90.in
# This module is used both by the CMake and the pre-existing build system
#
# copyright (c) 2017, Roberto Di Remigio.
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

import os
import sys
import subprocess
import re
import getpass
import platform
import datetime

def configure_file(rep, fname, **kwargs):
    ''' Configure a file.
    This functions acts more or less like CMake's configure_file command.

    :param rep:
       a {placeholder : replacement} dictionary
    :param fname:
       name of the file to be configured, without suffix
    :param \**kwargs:
       See below

    :Keyword arguments:
       * *in_path*  -- directory for the unconfigured file
       * *suffix*   -- suffix of the unconfigured file, with separators
       * *prefix*   -- prefix for the configured file
       * *out_path* -- directory for the configured file
    '''
    in_path = kwargs.get('in_path', os.getcwd())
    suffix = kwargs.get('suffix', '.in')
    out_path = kwargs.get('out_path', in_path)
    prefix = kwargs.get('prefix', '')
    fname_in = fname + suffix
    filedata = ''
    with open(os.path.join(in_path, fname_in), 'r') as fin:
        filedata = fin.read()
    rep = dict((re.escape(k), v) for k, v in rep.items())
    pattern = re.compile("|".join(list(rep.keys())))
    filedata = pattern.sub(lambda m: rep[re.escape(m.group(0))], filedata)
    fname_out = prefix + fname
    with open(os.path.join(out_path, fname_out), 'w+') as fout:
        fout.write(filedata)

def prepare_configuration_dictionary(**kwargs):
    '''
    Some of the configuration variables are set _via_ preprocessor variables.
    Others are set either in the configuration dictionary prepared by
    mkconfig.py or by variables set within the CMake build system.

    :Keyword arguments:
       * *build_type*
       * *Fortran_compiler*
       * *C_compiler*
       * *cmake_version*
       * *cmake_generator*
       * *mpi_launcher*
    '''

    build_type = kwargs.get('build_type', 'unknown')
    Fortran_compiler = kwargs.get('Fortran_compiler', 'unknown')
    C_compiler = kwargs.get('C_compiler', 'unknown')
    cmake_version = kwargs.get('cmake_version', 'Not built using CMake')
    cmake_generator = kwargs.get('cmake_generator', 'Not built using CMake')
    mpi_launcher = kwargs.get('mpi_launcher', 'unknown')

    conf_dict = {
        '@_user_name@' : getpass.getuser(),
        '@_host_name@' : platform.node(),
        '@_system@' : '{}-{}'.format(platform.system(), platform.release()),
        '@_cmake_version@' : cmake_version,
        '@_cmake_generator@' : cmake_generator,
        '@_build_type@' : build_type,
        '@_configuration_time@' : datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S [UTC]'),
        '@_python_version@' : '{}.{}.{}'.format(sys.version_info[0],
                                                sys.version_info[1],
                                                sys.version_info[2]),
        '@_Fortran_compiler@' : Fortran_compiler,
        '@_C_compiler@' : C_compiler,
        '@_pop_size@' : 'POP_SIZE',
        '@_det_size@' : 'DET_SIZE',
        '@_use_hdf5@' : 'DISABLE_HDF5',
        '@_use_lanczos@' : 'DISABLE_LANCZOS',
        '@_use_uuid@' : 'DISABLE_UUID',
        '@_use_scalapack@' : 'DISABLE_SCALAPACK',
        '@_use_single_precision@' : 'SINGLE_PRECISION',
        '@_use_popcnt@' : 'USE_POPCNT',
        '@_use_backtrace@' : 'DISABLE_BACKTRACE',
        '@_dsfmt_mexp@' : 'DSFMT_MEXP',
        '@_enable_mpi@' : 'PARALLEL',
        '@_mpi_launcher@' : mpi_launcher,
        '@_enable_omp@' : '_OPENMP',
        '@_git_last_commit_hash@' : 'unknown',
        '@_git_last_commit_author@' : 'unknown',
        '@_git_last_commit_date@' : 'unknown',
        '@_git_branch@' : 'unknown',
        '@_git_describe@' : '1.1-dev'
    }

    git_last_commit_hash = run_git('--no-pager show -s --pretty=format:%H -n 1')
    conf_dict.update({ '@_git_last_commit_hash@' : git_last_commit_hash })
    git_last_commit_author = run_git('--no-pager show -s --pretty=format:%aN -n 1')
    conf_dict.update({ '@_git_last_commit_author@' : git_last_commit_author })
    git_last_commit_date = run_git('--no-pager show -s --pretty=format:%ad -n 1')
    conf_dict.update({ '@_git_last_commit_date@' : git_last_commit_date })
    git_branch = run_git('rev-parse --abbrev-ref HEAD')
    conf_dict.update({ '@_git_branch@' : git_branch })
    git_describe = run_git('describe --abbrev=7 --long --always --dirty --tags')
    conf_dict.update({ '@_git_describe@' : git_describe })

    return conf_dict

def run_git(args):
    """Run Git with provided arguments"""
    try:
        p = subprocess.Popen('git ' + args,
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             shell=True)
        stdout_coded, stderr_coded = p.communicate()
        stdout = stdout_coded.decode('UTF-8')
        stderr = stderr_coded.decode('UTF-8')
        retcode = p.returncode
        if retcode < 0:
            sys.stderr.write('Git terminated by signal {}'.format(-retcode))
    except OSError as e:
        sys.stderr.write(stderr)
        sys.stderr.write('Git execution failed: {}'.format(e))
    return stdout.rstrip()
