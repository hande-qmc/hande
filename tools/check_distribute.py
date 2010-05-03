#!/usr/bin/python
'''Examine variables set by parsing input files and determine if they are being distributed to the other nodes or not.

This is exceptionally useful for sanity-checking when debugging parallel code.'''

import os
import re
import sys

# script resides in the tools subdirectory of the project.
tools_dir = os.path.dirname(sys.argv[0])
# hence location of the src subdirectory.
src_dir = os.path.abspath(os.path.join(tools_dir, '../src'))
# hence location of parse_input.F90
file = os.path.join(src_dir,'parse_input.F90')

f = open(file, 'r')

# Get keywords set in read_input.F90
start = re.compile('^ *do ! loop over lines in input file.', re.I)
end = re.compile('^ *end do ! end reading of input.', re.I)
read = re.compile('^.*? *call read[aif]|^.*? *call get[aif]', re.I)
setvar = re.compile('^ *[a-z_]+ ?=', re.I)
parentheses = re.compile('(?<=\()(.*?)(?=\(|\))')
data_in = False

variables = set([])

for line in f:
    if start.match(line):
        data_in = True
    if data_in:
        if read.match(line):
            # e.g. if (test) call readi(a)
            # obtain a.
            # 1. obtain 'readi(a)'
            fn_call = line.split('call')[-1].strip()
            # 2. obtain 'a'.
            var = parentheses.search(fn_call).group(0)
            variables.update([var])
        if setvar.match(line):
            # e.g. a = b 
            # obtain a.
            var = line.split('=')[-2].strip()
            variables.update([var])
    if end.match(line):
        data_in = False
        break

# special case: output filenames are not needed apart from on the head node.
variables.remove('hamiltonian_file')
variables.remove('determinant_file')

# Now get variables which are distributed in distribute_input
distributed = set([])
bcast = re.compile('(?<=call mpi_bcast\()(.*?)(?=,)', re.I)

for line in f:
    bcast_match = bcast.search(line)
    if bcast_match:
        distributed.update([bcast_match.group(0)])

# special case: option_set is used only for some book-keeping in distribute_input.
distributed.remove('option_set')

exit = 0
if distributed.difference(variables):
    print 'Distributed variables that are not read from input file are:', ' '.join(distributed.difference(variables))
    exit += 1
else:
    print 'All distributed variables can be set in the input file.'
if variables.difference(distributed):
    print 'Variables read from input file that are not distributed are:', ' '.join(variables.difference(distributed))
    exit += 2
else:
    print 'All variables set in the input file are distributed.'

sys.exit(exit)
