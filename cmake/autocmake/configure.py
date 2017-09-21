import os
import sys


def module_exists(module_name):
    try:
        __import__(module_name)
    except ImportError:
        return False
    else:
        return True


def check_cmake_exists(cmake_command):
    """
    Check whether CMake is installed. If not, print
    informative error message and quits.
    """
    from subprocess import Popen, PIPE

    p = Popen('{0} --version'.format(cmake_command),
              shell=True,
              stdin=PIPE,
              stdout=PIPE)
    if not ('cmake version' in p.communicate()[0].decode('UTF-8')):
        sys.stderr.write('   This code is built using CMake\n\n')
        sys.stderr.write('   CMake is not found\n')
        sys.stderr.write('   get CMake at http://www.cmake.org/\n')
        sys.stderr.write('   on many clusters CMake is installed\n')
        sys.stderr.write('   but you have to load it first:\n')
        sys.stderr.write('   $ module load cmake\n')
        sys.exit(1)


def setup_build_path(build_path):
    """
    Create build directory. If this already exists, print informative
    error message and quit.
    """
    if os.path.isdir(build_path):
        fname = os.path.join(build_path, 'CMakeCache.txt')
        if os.path.exists(fname):
            sys.stderr.write('aborting setup\n')
            sys.stderr.write('build directory {0} which contains CMakeCache.txt already exists\n'.format(build_path))
            sys.stderr.write('remove the build directory and then rerun setup\n')
            sys.exit(1)
    else:
        os.makedirs(build_path, 0o755)


def test_adapt_cmake_command_to_platform():

    cmake_command = "FC=foo CC=bar CXX=RABOOF cmake -DTHIS -DTHAT='this and that cmake' .."
    res = adapt_cmake_command_to_platform(cmake_command, 'linux')
    assert res == cmake_command
    res = adapt_cmake_command_to_platform(cmake_command, 'win32')
    assert res == "set FC=foo && set CC=bar && set CXX=RABOOF && cmake -DTHIS -DTHAT='this and that cmake' .."

    cmake_command = "cmake -DTHIS -DTHAT='this and that cmake' .."
    res = adapt_cmake_command_to_platform(cmake_command, 'linux')
    assert res == cmake_command
    res = adapt_cmake_command_to_platform(cmake_command, 'win32')
    assert res == cmake_command


def adapt_cmake_command_to_platform(cmake_command, platform):
    """
    Adapt CMake command to MS Windows platform.
    """
    if platform == 'win32':
        pos = cmake_command.find('cmake')
        s = ['set {0} &&'.format(e) for e in cmake_command[:pos].split()]
        s.append(cmake_command[pos:])
        return ' '.join(s)
    else:
        return cmake_command


def run_cmake(command, build_path, default_build_path):
    """
    Execute CMake command.
    """
    from subprocess import Popen, PIPE
    from shutil import rmtree

    topdir = os.getcwd()
    os.chdir(build_path)
    p = Popen(command,
              shell=True,
              stdin=PIPE,
              stdout=PIPE,
              stderr=PIPE)
    stdout_coded, stderr_coded = p.communicate()
    stdout = stdout_coded.decode('UTF-8')
    stderr = stderr_coded.decode('UTF-8')

    # print cmake output to screen
    print(stdout)

    if stderr:
        # we write out stderr but we do not stop yet
        # this is because CMake warnings are sent to stderr
        # and they might be benign
        sys.stderr.write(stderr)

    # write cmake output to file
    with open('cmake_output', 'w') as f:
        f.write(stdout)

    # change directory and return
    os.chdir(topdir)

    # to figure out whether configuration was a success
    # we check for 3 sentences that should be part of stdout
    configuring_done = '-- Configuring done' in stdout
    generating_done = '-- Generating done' in stdout
    build_files_written = '-- Build files have been written to' in stdout
    configuration_successful = configuring_done and generating_done and build_files_written

    if configuration_successful:
        save_setup_command(sys.argv, build_path)
        print_build_help(build_path, default_build_path)
    else:
        if (build_path == default_build_path):
            # remove build_path iff not set by the user
            # otherwise removal can be dangerous
            rmtree(default_build_path)


def print_build_help(build_path, default_build_path):
    """
    Print help text after configuration step is done.
    """
    print('   configure step is done')
    print('   now you need to compile the sources:')
    if (build_path == default_build_path):
        print('   $ cd build')
    else:
        print('   $ cd ' + build_path)
    print('   $ make')


def save_setup_command(argv, build_path):
    """
    Save setup command to a file.
    """
    file_name = os.path.join(build_path, 'setup_command')
    with open(file_name, 'w') as f:
        f.write(' '.join(argv[:]) + '\n')


def configure(root_directory, build_path, cmake_command, only_show):
    """
    Main configure function.
    """
    default_build_path = os.path.join(root_directory, 'build')

    # check that CMake is available, if not stop
    check_cmake_exists('cmake')

    # deal with build path
    if build_path is None:
        build_path = default_build_path
    if not only_show:
        setup_build_path(build_path)

    cmake_command = adapt_cmake_command_to_platform(cmake_command, sys.platform)

    print('{0}\n'.format(cmake_command))
    if only_show:
        sys.exit(0)

    run_cmake(cmake_command, build_path, default_build_path)
