#!/usr/bin/env python2
'''Usage: extract.py filename

Extract data from the output of a hubbard_fciqmc calculation.'''

import optparse
import re
import sys

# If using python <2.5, any doesn't exist. Define our own.
try:
    any([])
except NameError:
    def any(iterable):
        '''any(iterable) -> bool
        
Return True if bool(x) is True for any x in the iterable.
'''
        for val in iterable:
            if val:
                return True
        return False

# data elements to find.
class DataElement(object):
    '''Holder for a data element which is searched for in an output file.'''
    def __init__(self, name, regex, line_index, regex_flags=None, verbose=False):
        '''Initialise a DataElement object.

name: name of data item.  Used in output table.  Should not contain spaces for testcode compatibility.
regex: regular expression used to find the data element.  re.match is used, so the regular expression should start from the beginning of the line.
line_index: the index of the data value in the line, when split by spaces.
regex_flags (optional): flags to use in the compiled regular expression.
verbose (optional): flag which if true denotes data should only be included in verbose output.  Dafeult: False.
'''
        self.name = name
        if regex_flags:
            self.regex = re.compile(regex, regex_flags)
        else:
            self.regex = re.compile(regex)
        self.line_index = line_index
        self.value = None
        self.verbose = verbose
    def fmt(self, padding=0):
        '''Return a format string which will hold both the name and value of the data item without truncation.

padding (optional integer): amount of space to add to format string.
'''
        if self.value:
            return '%%-%is' % (max(len(str(self.value)), len(str(self.name)))+padding)
        else:
            return '%%-%is' % (len(str(self.name))+padding)

ALL_DATA_ELEMENTS = [
        DataElement('exact_ground_state','^ *Exact ground state:', -1),
        DataElement('lanczos_ground_state','^ *Lanczos ground state:', -1),
        DataElement('E0','^ *E0 = <D0|H|D0> =', -1),
        DataElement('final_shift','^ final shift =', -1),
        DataElement('final_proj_E','^ final proj. energy =', -1),
        DataElement('E0+final_shift','^ E0 \+ shift =', -1, verbose=True),
        DataElement('hilbert_space_size','^ Monte-Carlo estimate of size of space is:', -1),
    ]

def parse_options(args):
    '''Parse command line options.'''

    parser = optparse.OptionParser(usage = __doc__)
    parser.add_option('-v', '--verbose', action='store_true', default=False, help='Show all extracted data (including data items which are sums of other data items. Default: %default.')
    (options, filename) = parser.parse_args(args)

    if len(filename) != 1:
        parser.print_help()
        sys.exit(1)

    return options, filename[0]

def extract(data_elements, filename):
    '''Extract the data from the file specified.

data_elements: list of DataElement objects.
filename: file to examine.
'''

    f = open(filename)

    for line in f:
        for data_item in data_elements:
            if re.match(data_item.regex, line):
                data_item.value = line.split()[data_item.line_index]

    f.close()

    padding = 3

    if any(d.value for d in data_elements): 
        # Have data to output.

        # print header
        for d in data_elements:
            if d.value:
                print d.fmt(padding) % d.name,
        print

        # print values
        for d in data_elements:
            if d.value:
                print d.fmt(padding) % d.value,
        print

    return None

if __name__ == '__main__':
    (options, filename) = parse_options(sys.argv[1:])
    if options.verbose:
        data_elements = ALL_DATA_ELEMENTS
    else:
        data_elements = [d for d in ALL_DATA_ELEMENTS if not d.verbose]
    extract(data_elements, filename)
