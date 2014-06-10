'''Data extraction from the (standard) output of a HANDE calculation.'''

import numpy
import pandas as pd
import os
import re
import sys
import tempfile

def extract_data_sets(filenames, temp_file=True):
    '''Extract QMC data tables from multiple HANDE calculations.

Parameters
----------
filenames : list of strings
    names of files containing HANDE QMC calculation output.
temp_file : bool
    passed to extract_data.

Returns
-------
metadata : :class:`pandas.DataFrame`
    metadata (i.e. calculation information, parameters and settings) extracted
    from output files.
data : :class:`pandas.DataFrame`
    HANDE QMC data.  Each calculation (and each beta loop in the case of DMQMC)
    is labelled separately using a hierarchical index.

See Also
--------
``extract_data`` : underlying data extraction implementation.
'''

    offset = 0
    data = []
    metadata = []
    for filename in filenames:
        (calc_metadata, calc_data) = extract_data(filename, offset, temp_file)
        data.append(calc_data)
        metadata.append(calc_metadata)
        offset += data[-1].index.levshape[0]
    return (pd.concat(metadata).unstack(level=0), pd.concat(data))

def extract_data(filename, offset=0, temp_file=True):
    '''Extract QMC data table from a HANDE calculation.

.. note::

    We assume that the data table is continuous (i.e. not split by anything
    other than comments) and that the final line of the table is within the
    final 2KB of the file.

Parameters
----------
filename : string
    name of file containing the HANDE QMC calculation output.
offset : int
    first index for the calculation level label in data (default: 0).
temp_file : bool
    if True convert data tables in large files (ie those over 8MB) to CSV format
    before parsing.  This is much faster and has lower memory usage as it allows
    pandas' C parser to be used.  This requires a temporary file to be created
    in TMPDIR (default: /tmp).

Returns
-------
metadata : :class:`pandas.Series`
    metadata (i.e. calculation information, parameters and settings) extracted
    from output file..
data : :class:`pandas.DataFrame`
    HANDE QMC data.  Multiple loops over iterations (e.g. in DMQMC) are labelled
    using a hierarchical index, starting from the supplied offset, otherwise the
    entire data set is given the offset as the hierarchical index.
'''

    md_regex = dict(
        UUID = '^ Calculation UUID:',
        calc = '^ *(fciqmc|ccmc|ifciqmc|iccmc|dmqmc|idmqmc) *$',
        sym = r'\bsym\b +\d+',
        ms = 'ms +-*\d+',
        nel = 'nel',
        nbasis = 'Number of basis functions:',
        truncation = 'truncation_level',
        tau = 'tau',
        seed = 'random number generator with a seed of',
        bit_length = 'Bit-length',
        ref = 'Reference determinant, \|D0> = ',
        ref_energy = r'E0 = <D0\|H\|D0>',
        psingle = 'Probability of attempting a single',
        pdouble = 'Probability of attempting a double',
        init_pop = 'Initial population on',
        git_hash = 'VCS BASE repository version:',
        target = 'varyshift_target',
        shift_damping = 'shift_damping',
        mc_cycles = 'mc_cycles'
    )
    md_int = 'sym ms nel nbasis truncation seed bit_length'.split()
    md_float = 'tau ref_energy psingle pdouble init_pop'.split()
    for (k,v) in md_regex.items():
        md_regex[k] = re.compile(v, re.IGNORECASE)

    # Read metadata and figure out how the start line of the data table.
    f = open(filename, 'r')
    start_line = 0
    metadata = pd.Series(index=md_regex.keys(), name='metadata')
    unseen_calc = True
    have_git_hash_next = False
    for line in f:
        start_line += 1
        # extract metadata.
        if have_git_hash_next:
             metadata['git_hash'] = line.split()[0]
             have_git_hash_next = False
        for (k,v) in md_regex.items():
            if v.search(line):
                # Special cases for unusual formats...
                if k == 'ref':
                    metadata[k] = v.split(line)[-1].strip()
                elif k == 'calc' and unseen_calc:
                    unseen_calc = False
                    metadata[k] = line.split()[0]
                elif k == 'seed':
                    metadata[k] = int(float(line.split()[-1]))
                elif k == 'git_hash':
                    have_git_hash_next = True
                else:
                    val = line.split()[-1]
                    if k in md_int:
                        metadata[k] = int(val)
                    elif k in md_float:
                        metadata[k] = float(val)
                    else:
                        metadata[k] = val
        # Hunt for start of data table.
        if ' # iterations' in line:
            # Columns are separated by at least two spaces but each column name
            # can contain words separated by just one space.
            column_names = re.split('   *', line[3:].strip())
            break
    f.close()

    # Hunt for the end of the table.
    skip_footer = 0
    end_lines = _get_last_lines(filename)
    for line in end_lines[::-1]:
        if re.match('  *\d\d?', line):
            break
        else:
            skip_footer += 1

    if float(os.path.getsize(filename))/1024 < 8000 or not temp_file:
        # Read table --- only read the first N columns, where N is the number of
        # column names found.
        data = pd.io.parsers.read_table(filename, sep='\s+', engine='python',
                   skiprows=start_line, skipfooter=skip_footer,
                   names=column_names, comment='#')
        # Remove comment lines and convert all columns to numeric data.
        # Lines starting with a comment have been set to NaN in the iterations
        # column.
        try:
            data.dropna(subset=['iterations'], inplace=True)
        except TypeError:
            # Be slightly less efficient if using pandas version < 0.13.
            data = data.dropna(subset=['iterations'])
        data.reset_index(drop=True, inplace=True)
    else:
        # Work around pandas slow and very memory-hungry pure-python parser by
        # converting the data table into CSV format (and stripping out comments whilst we're at it)
        # and reading that in.
        tmp_csv = _convert_to_csv(filename, start_line)
        data = pd.io.parsers.read_csv(tmp_csv, names=column_names)
        os.remove(tmp_csv)

    data = data.convert_objects(convert_numeric=True, copy=False)

    unique_iterations = data['iterations'].unique()
    # if repeated iterations:
    #     assume there's a full calculation (i.e. DMQMC style) followed by
    #     a subsequent ones (the final one might not be complete), e.g.
    #         iteration  label_1 label_2 ...
    #            10       ....
    #            20       ....
    #            ..       ....
    #           300       ....
    #            10       ....
    #            20       ....
    #            ..       ....
    #           300       ....
    #            10       ....
    #            20       ....
    #     then label each set of iterations with X.
    # Otherwise just label the data just read in by the same X.
    niterations = len(data)
    nunique = len(unique_iterations)
    nitems = int(numpy.ceil(float(niterations)/nunique))
    data = [ data[i*nunique:(i+1)*nunique] for i in range(nitems) ]
    multi_keys = [i+offset for i in range(len(data))]
    data = pd.concat(data, keys=multi_keys)
    data.index.names = ['calc', '']
    metadata = pd.concat([metadata], keys=[multi_keys[0]])
    metadata.index.names = ['calc', '']

    # Do we have an old table?  If so, rename the headings to the new ones for
    # convenience...
    data.rename(inplace=True, columns={
            'Instant shift': 'Shift',
            '\sum H_0j Nj': '\sum H_0j N_j',
            '# D0': 'N_0',
            '# particles': '# H psips',
        })

    return (metadata, data)

def _get_last_lines(filename, bytes=2048):
    '''Get the lines within a given number of bytes from the end of the file.

Parameters
----------
filename : string
    name of file.
bytes : int
    number of bytes from the end of file to examine.
    
Returns
-------
last_lines : list
    list of lines.
'''

    f = open(filename, 'r')
    f.seek (0, 2)
    fsize = f.tell()
    f.seek(max(fsize-2048, 0), 0) # Get the last 2k of file.
    lines = f.readlines()
    f.close()
    return lines

def _convert_to_csv(filename, start_line, comment='#'):
    '''Convert the HANDE data table in a file to CSV format.

Parameters
----------
filename : string
    name of HANDE output file.
start_line : int
    line number at which the QMC data table starts.
comment : string
    Single character which indicates a comment line if its the first
    non-whitespace character in a line.

.. note::

    We assume the table finishes at the first blank line after start_line.

Returns
-------
temp_file_name : string
    name of the CSV temporary file.
'''
    tmp = tempfile.NamedTemporaryFile(delete=False, mode='w')
    nlines = 0
    f = open(filename)
    for line in f:
        nlines += 1
        if nlines >= start_line:
            csv_line = ','.join(line.strip().split())
            if not csv_line:
                # blank line => end of data table.
                break
            elif not csv_line.startswith(comment):
                tmp.write(csv_line+'\n')
    f.close()
    tmp.close()

    return tmp.name
