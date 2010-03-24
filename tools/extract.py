#!/usr/bin/env python
'''Extract data from the output of a hubbard_fciqmc calculation.

Usage: extract.py filename
'''

import re
import sys

class DataElement(object):
    '''Holder for a data element which is searched for in an output file.'''
    def __init__(self, name, regex, line_index, regex_flags=None):
        '''Initialise a DataElement object.

name: name of data item.  Used in output table.  Should not contain spaces for testcode compatibility.
regex: regular expression used to find the data element.  re.match is used, so the regular expression should start from the beginning of the line.
line_index: the index of the data value in the line, when split by spaces.
regex_flags (optional): flags to use in the compiled regular expression.
'''
        self.name = name
        if regex_flags:
            self.regex = re.compile(regex, regex_flags)
        else:
            self.regex = re.compile(regex)
        self.line_index = line_index
        self.value = None
    def fmt(self, padding=0):
        '''Return a format string which will hold both the name and value of the data item without truncation.

padding (optional integer): amount of space to add to format string.
'''
        if self.value:
            return '%%-%is' % (max(len(str(self.value)), len(str(self.name)))+padding)
        else:
            return '%%-%is' % (len(str(self.name))+padding)

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
data_elements = [
        DataElement('exact_ground_state','^ *Exact ground state:', -1),
        DataElement('lanczos_ground_state','^ *Lanczos ground state:', -1),
        DataElement('E0','^ *E0 = <D0|H|D0> =', -1),
        DataElement('final_shift','^ final shift =', -1),
        DataElement('av_shift','^ av\. shift =', -1),
        DataElement('final_proj_E','^ final proj. energy =', -1),
        DataElement('av_proj_E','^ av\. proj. energy =', -1),
        DataElement('E0+final_shift','^ E0 \+ shift =', -1),
        DataElement('E0+av_shift','^ E0 \+ av\. shift =', -1),
        DataElement('E0+final_proj_E','^ E0 \+ av\. proj\. energy =', -1),
    ]

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
    if len(sys.argv) == 2:
        filename = sys.argv[1]
        extract(data_elements, filename)
    else:
        print __doc__,
