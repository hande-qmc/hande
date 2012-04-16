#!/usr/bin/env python
'''A simple scipt to check that all the input options are tested and documented.'''

import glob
import os
import re
import sys

# script resides in the tools subdirectory of the project.
TOOLS_DIR = os.path.dirname(sys.argv[0])
# hence location of the src subdirectory.
SRC_DIR = os.path.abspath(os.path.join(TOOLS_DIR, '../src'))
# hence location of the documentation subdirectory.
DOCS_DIR = os.path.abspath(os.path.join(TOOLS_DIR, '../documentation/users_guide'))
# hence location of the test suite subdirectory.
TEST_DIR = os.path.abspath(os.path.join(TOOLS_DIR, '../test_suite'))

# glob for documentation files
DOC_FILES = os.path.join(DOCS_DIR, '*.rst')
# input parser source file
INPUT_PARSER = os.path.join(SRC_DIR, 'parse_input.F90')

# glob for test input files
TEST_FILES = os.path.join(TEST_DIR, '*/*.inp')

# useful, common, regexes.
COMMENT = re.compile('^ *!')

def get_keywords(filename):
    '''Get keywords used by input parser, which are given in the form

case('KEYWORD')

or

case('KEYWORD1','KEYWORD2')

where the latter two keywords are synonyms.

Return a list of all keywords and a list of lists of equivalent keywords.'''

    all_keywords = []
    unique_keywords = []

    f = open(filename)
    keywords = re.compile('^ *case\(".*"\)|^ *case\(\'.*\'\)', re.I)
    strip = re.compile('!.*$|case\(|\'|"|\)|  *|\n')
    for line in f:
        if keywords.match(line):
            # remove comments
            line = re.sub('!.*$','',line)
            # remove case, brackets, quote marks and whitespace
            line = strip.sub('', line)
            words = line.split(',')
            all_keywords.extend(words)
            unique_keywords.append(words)

    return (all_keywords,unique_keywords)

def check_docs(doc_files, keywords):
    '''Search for keywords in doc_files and print those not found.'''

    doc_files = glob.glob(doc_files)

    for f in doc_files:
        if not os.path.isdir(os.path.dirname(f)):
            print('Documentation directory does not exist.')
    if not doc_files:
        print('Cannot find documentation files.')

    keywords = [(k,re.compile(r'\*\*%s\*\*' % k,re.I)) for k in keywords]
    for filename in doc_files:
        not_found = []
        f = open(filename)
        docs = ''.join(f.readlines())
        for k in keywords:
            if not k[1].search(docs):
                not_found.append(k)
        keywords = not_found
        f.close()

    if not_found:
        print('Undocumented keywords are:')
        print()
        for k in keywords:
            print(k[0])
        print()
    else:
        print('All keywords appear in the documentation.')

def check_tests(test_files, keywords):
    '''Search for keywords in test_files and print those not found.'''

    test_files = glob.glob(test_files)

    for f in test_files:
        if not os.path.isdir(os.path.dirname(f)):
            print('Test suite directory does not exist.')
    if not test_files:
        print('Cannot find test suite input files.')

    patterns = []
    for k in keywords:
        pattern = '|'.join([r'\b%s\b' % word for word in k])
        patterns.append((k,re.compile(pattern,re.I)))
    keywords = patterns
    for filename in test_files:
        not_found = []
        f = open(filename)
        inp = ' '.join(f.readlines())
        for k in keywords:
            if not k[1].search(inp):
                not_found.append(k)
        keywords = not_found
        f.close()

    if not_found:
        print('Untested keywords are:')
        print()
        for k in not_found:
            print(' or '.join(k[0]))
        print()
    else:
        print('All keywords (or their synonyms) appear in the input files in the test suite.')

if __name__ == '__main__':

    (all_keywords,unique_keywords) = get_keywords(INPUT_PARSER)

    check_docs(DOC_FILES, all_keywords)

    check_tests(TEST_FILES, unique_keywords)
