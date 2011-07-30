#!/usr/bin/env python2
'''Examine variables set by parsing input files and determine if they are being distributed to the other nodes or not.

This is exceptionally useful for sanity-checking when debugging parallel code.'''

import os
import re
import sys

# script resides in the tools subdirectory of the project.
TOOLS_DIR = os.path.dirname(sys.argv[0])
# hence location of the src subdirectory.
SRC_DIR = os.path.abspath(os.path.join(TOOLS_DIR, '../src'))

# useful, common, regexes.
BCAST = re.compile('(?<=call mpi_bcast\()(.*?)(?=,)', re.I)
SCATTERV = re.compile('(?<=call mpi_scatterv\()(.*?)(?=,)', re.I)
COMMENT = re.compile('^ *!')

def check_distributed(filename, variables, distributed):

    exit = 0
    if distributed.difference(variables):
        print 'Distributed variables that are not read in file %s are: %s.' % (filename, ', '.join(distributed.difference(variables)))
        exit += 1
    else:
        print 'All distributed variables can be set in the file %s.' % (filename)
    if variables.difference(distributed):
        print 'Variables read in file %s that are not distributed are: %s.' % (filename, ', '.join(variables.difference(distributed)))
        exit += 2
    else:
        print 'All variables set in the input file are distributed.'
    return exit

def test_input(input_file):

    # hence location of input_file
    file = os.path.join(SRC_DIR, input_file)

    f = open(file, 'r')

    # Get keywords set in source file.
    start = re.compile('^ *do ! loop over lines', re.I)
    end = re.compile('^ *end do ! end reading.', re.I)
    read = re.compile('^.*? *call read[aif]|^.*? *call get[aif]', re.I)
    setvar = re.compile('^ *[a-z_]+ ?=', re.I)
    parentheses = re.compile('(?<=\()(.*?)(?=\(|\))')
    data_in = False

    variables = set([])

    for line in f:
        if start.match(line):
            data_in = True
        if data_in:
            if not COMMENT.match(line):
                if read.match(line):
                    # e.g. if (test) call readi(a)
                    # obtain a.
                    # 1. obtain 'readi(a)'
                    fn_call = line.split('call')[-1].strip()
                    # 2. obtain 'a'.
                    var = parentheses.search(fn_call).group(0)
                    variables.update([var])
                elif setvar.match(line):
                    # e.g. a = b 
                    # obtain a.
                    var = line.split('=')[-2].strip()
                    variables.update([var])
        if end.match(line):
            data_in = False
            break

    # Now get variables which are distributed in the source file.
    distributed = set([])

    for line in f:
        bcast_match = BCAST.search(line)
        if bcast_match and not COMMENT.match(line):
            distributed.update([bcast_match.group(0)])

    # special case: output filenames are not needed apart from on the head node, reading the restart file is only read on the head node.
    for v in ['hamiltonian_file', 'determinant_file', 'binary_fmt_in', 'binary_fmt_out']:
        if v in variables:
            variables.remove(v)
    # special case: option_set is used only for some book-keeping in distribute_input; comms_read is similarly used in fciqmc_interact.
    for v in ['option_set', 'comms_read']:
        if v in distributed:
            distributed.remove(v)

    return check_distributed(input_file, variables, distributed)

def test_restart():

    # hence location of restart.F90
    input_file = 'restart.F90'
    file = os.path.join(SRC_DIR, input_file)

    f = open(file, 'r')

    # Get keywords set in restart.F90
    start = re.compile('^ *subroutine read_restart', re.I)
    end = re.compile('^ *end subroutine read_restart', re.I)
    read = re.compile('(?<=call read_in\(io,)(.*?)(?=\))', re.I)
    # remove commas and array slices
    data_extract = re.compile(',|(\(.*?\))')
    data_in = False

    variables = set([])
    distributed = set([])

    # Get values read in and distributed.
    for line in f:
        if start.match(line):
            data_in = True
        if data_in:
            if not COMMENT.match(line):
                read_search = read.search(line)
                if read_search:
                    data = read_search.group(0).strip().split(',')
                    data = [re.sub(data_extract, '', d) for d in data]
                    variables.update(data)
                else:
                    for mpi_regex in [BCAST, SCATTERV]:
                        mpi_match = mpi_regex.search(line)
                        if mpi_match:
                            distributed.update([mpi_match.group(0)])
        if end.match(line):
            data_in = False

    # Remove special cases...
    for v in ['det', 'pop', 'energy', 'junk', 'tot_walkers']:
        variables.remove(v)
    for v in ['done', 'scratch_energies', 'spawned_walkers']:
        distributed.remove(v)

    return check_distributed(input_file, variables, distributed)

exit = test_input('parse_input.F90') + 10*test_input('interact.F90') + 100*test_restart()

sys.exit(exit)
