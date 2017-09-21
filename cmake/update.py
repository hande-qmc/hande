#!/usr/bin/env python

import os
import sys


if sys.version_info[0] == 2 and sys.version_info[1] < 7:
    sys.stderr.write("ERROR: update.py requires at least Python 2.7\n")
    sys.exit(-1)


AUTOCMAKE_GITHUB_URL = 'https://github.com/coderefinery/autocmake/raw/master/'


def check_for_yaml():
    try:
        import yaml
    except:
        sys.stderr.write("ERROR: you need to install the pyyaml package\n")
        sys.exit(-1)


def print_progress_bar(text, done, total, width):
    """
    Print progress bar.
    """
    if total > 0:
        n = int(float(width) * float(done) / float(total))
        sys.stdout.write("\r{0} [{1}{2}] ({3}/{4})".format(text, '#' * n, ' ' * (width - n), done, total))
        sys.stdout.flush()


def flat_add(l, x):
    if isinstance(x, int):
        l.append(x)
        return l
    elif isinstance(x, str):
        l.append(x)
        return l
    else:
        return l + x


def fetch_modules(config, relative_path, download_directory):
    """
    Assemble modules which will
    be included in CMakeLists.txt.
    """
    from collections import Iterable, namedtuple, defaultdict
    from autocmake.extract import extract_list, to_d, to_l
    from autocmake.parse_rst import parse_cmake_module

    cleaned_config = defaultdict(lambda: [])

    modules = []
    Module = namedtuple('Module', 'path name')

    num_sources = len(extract_list(config, 'source'))

    print_progress_bar(text='- assembling modules:',
                       done=0,
                       total=num_sources,
                       width=30)

    if 'modules' in config:
        i = 0
        for t in config['modules']:
            for k, v in t.items():

                d = to_d(v)
                for _k, _v in to_d(v).items():
                    cleaned_config[_k] = flat_add(cleaned_config[_k], _v)

                # fetch sources and parse them
                if 'source' in d:
                    for src in to_l(d['source']):
                        i += 1

                        # we download the file
                        module_name = os.path.basename(src)
                        if 'http' in src:
                            path = download_directory
                            name = 'autocmake_{0}'.format(module_name)
                            dst = os.path.join(download_directory, 'autocmake_{0}'.format(module_name))
                            fetch_url(src, dst)
                            file_name = dst
                            fetch_dst_directory = download_directory
                        else:
                            if os.path.exists(src):
                                path = os.path.dirname(src)
                                name = module_name
                                file_name = src
                                fetch_dst_directory = path
                            else:
                                sys.stderr.write("ERROR: {0} does not exist\n".format(src))
                                sys.exit(-1)

                        # we infer config from the module documentation
                        # dictionary d overrides the configuration in the module documentation
                        # this allows to override interpolation inside the module
                        with open(file_name, 'r') as f:
                            parsed_config = parse_cmake_module(f.read(), d)
                            for _k2, _v2 in parsed_config.items():
                                if _k2 not in to_d(v):
                                    # we add to clean_config only if the entry does not exist
                                    # in parent autocmake.yml already
                                    # this allows to override
                                    cleaned_config[_k2] = flat_add(cleaned_config[_k2], _v2)

                        modules.append(Module(path=path, name=name))
                        print_progress_bar(text='- assembling modules:',
                                           done=i,
                                           total=num_sources,
                                           width=30)
        print('')

    return modules, cleaned_config


def process_yaml(argv):
    from autocmake.parse_yaml import parse_yaml
    from autocmake.generate import gen_cmakelists, gen_setup
    from autocmake.extract import extract_list

    project_root = argv[1]
    if not os.path.isdir(project_root):
        sys.stderr.write("ERROR: {0} is not a directory\n".format(project_root))
        sys.exit(-1)

    # read config file
    print('- parsing autocmake.yml')
    with open('autocmake.yml', 'r') as stream:
        config = parse_yaml(stream)

    if 'name' in config:
        project_name = config['name']
    else:
        sys.stderr.write("ERROR: you have to specify the project name in autocmake.yml\n")
        sys.exit(-1)
    if ' ' in project_name.rstrip():
        sys.stderr.write("ERROR: project name contains a space\n")
        sys.exit(-1)

    if 'language' in config:
        project_language = ' '.join(config['language']) if isinstance(config['language'], list) else config['language']
    else:
        sys.stderr.write("ERROR: you have to specify the project language(s) in autocmake.yml\n\n")
        sys.stderr.write("# for instance like this (several languages):\nlanguage:\n  - CXX\n  - Fortran\n\n")
        sys.stderr.write("# or like this (one language):\nlanguage: Fortran\n\n")
        sys.exit(-1)

    if 'min_cmake_version' in config:
        min_cmake_version = config['min_cmake_version']
    else:
        sys.stderr.write("ERROR: you have to specify min_cmake_version in autocmake.yml\n")
        sys.exit(-1)

    if 'default_build_type' in config:
        default_build_type = config['default_build_type'].lower()
    else:
        sys.stderr.write("ERROR: you have to specify default_build_type in autocmake.yml\n\n")
        sys.stderr.write("# for instance like this (debug, release, relwithdebinfo, or minsizerel):\ndefault_build_type: release\n\n")
        sys.exit(-1)

    if 'setup_script' in config:
        setup_script_name = config['setup_script']
    else:
        setup_script_name = 'setup'

    # get relative path from setup script to this directory
    relative_path = os.path.relpath(os.path.abspath('.'), project_root)

    download_directory = 'downloaded'
    if not os.path.exists(download_directory):
        os.makedirs(download_directory)

    # fetch modules from the web or from relative paths
    modules, cleaned_config = fetch_modules(config, relative_path, download_directory)

    # fetch files which are not parsed
    for src in cleaned_config['fetch']:
        dst = os.path.join(download_directory, os.path.basename(src))
        fetch_url(src, dst)

    # print warnings
    for warning in cleaned_config['warning']:
        print('- WARNING: {0}'.format(warning))

    # create CMakeLists.txt
    print('- generating CMakeLists.txt')
    s = gen_cmakelists(project_name, project_language, min_cmake_version, default_build_type, relative_path, modules)
    with open(os.path.join(project_root, 'CMakeLists.txt'), 'w') as f:
        f.write('{0}\n'.format('\n'.join(s)))

    # create setup script unless it is 'None' or 'none'
    if setup_script_name.lower() != 'none':
        print('- generating setup script')
        s = gen_setup(cleaned_config, default_build_type, relative_path, setup_script_name)
        file_path = os.path.join(project_root, setup_script_name)
        with open(file_path, 'w') as f:
            f.write('{0}\n'.format('\n'.join(s)))
        if sys.platform != 'win32':
            make_executable(file_path)


def main(argv):
    """
    Main function.
    """

    if len(argv) != 2:
        sys.stderr.write("\nYou can update a project in two steps.\n\n")
        sys.stderr.write("Step 1: Update or create infrastructure files\n")
        sys.stderr.write("        which will be needed to configure and build the project:\n")
        sys.stderr.write("        $ {0} --self\n\n".format(argv[0]))
        sys.stderr.write("Step 2: Create CMakeLists.txt and setup script in PROJECT_ROOT:\n")
        sys.stderr.write("        $ {0} <PROJECT_ROOT>\n".format(argv[0]))
        sys.stderr.write("        example:\n")
        sys.stderr.write("        $ {0} ..\n".format(argv[0]))
        sys.exit(-1)

    if argv[1] in ['-h', '--help']:
        print('Usage:')
        for t, h in [('python update.py --self',
                      'Update this script and fetch or update infrastructure files under autocmake/.'),
                     ('python update.py <builddir>',
                      '(Re)generate CMakeLists.txt and setup script and fetch or update CMake modules.'),
                     ('python update.py (-h | --help)',
                      'Show this help text.')]:
            print('  {0:30} {1}'.format(t, h))
        sys.exit(0)

    if argv[1] == '--self':
        # update self
        if not os.path.isfile('autocmake.yml'):
            print('- fetching example autocmake.yml')
            fetch_url(
                src='{0}example/autocmake.yml'.format(AUTOCMAKE_GITHUB_URL),
                dst='autocmake.yml'
            )
        if not os.path.isfile('.gitignore'):
            print('- creating .gitignore')
            with open('.gitignore', 'w') as f:
                f.write('*.pyc\n')
        for f in ['autocmake/configure.py',
                  'autocmake/__init__.py',
                  'autocmake/external/docopt.py',
                  'autocmake/external/__init__.py',
                  'autocmake/generate.py',
                  'autocmake/extract.py',
                  'autocmake/interpolate.py',
                  'autocmake/parse_rst.py',
                  'autocmake/parse_yaml.py',
                  'update.py']:
            print('- fetching {0}'.format(f))
            fetch_url(
                src='{0}{1}'.format(AUTOCMAKE_GITHUB_URL, f),
                dst='{0}'.format(f)
            )
        sys.exit(0)

    process_yaml(argv)


def make_executable(path):
    # http://stackoverflow.com/a/30463972
    mode = os.stat(path).st_mode
    mode |= (mode & 0o444) >> 2    # copy R bits to X
    os.chmod(path, mode)


def fetch_url(src, dst):
    """
    Fetch file from URL src and save it to dst.
    """
    # we do not use the nicer sys.version_info.major
    # for compatibility with Python < 2.7
    if sys.version_info[0] > 2:
        import urllib.request

        class URLopener(urllib.request.FancyURLopener):
            def http_error_default(self, url, fp, errcode, errmsg, headers):
                sys.stderr.write("ERROR: could not fetch {0}\n".format(url))
                sys.exit(-1)
    else:
        import urllib

        class URLopener(urllib.FancyURLopener):
            def http_error_default(self, url, fp, errcode, errmsg, headers):
                sys.stderr.write("ERROR: could not fetch {0}\n".format(url))
                sys.exit(-1)

    dirname = os.path.dirname(dst)
    if dirname != '':
        if not os.path.isdir(dirname):
            os.makedirs(dirname)

    opener = URLopener()
    opener.retrieve(src, dst)


if __name__ == '__main__':
    check_for_yaml()
    main(sys.argv)
