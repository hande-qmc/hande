'''Data extraction from the (standard) output of a HANDE calculation.'''

import numpy
import pandas as pd
import re
import sys

def extract_data_sets(filenames):
    '''Extract QMC data tables from multiple HANDE calculations.

Parameters
----------
filenames : list of strings
    names of files containing HANDE QMC calculation output.

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
        (calc_metadata, calc_data) = extract_data(filename, offset)
        data.append(calc_data)
        metadata.append(calc_metadata)
        offset += data[-1].index.levshape[0]
    return (pd.concat(metadata).unstack(level=0), pd.concat(data))

def extract_data(filename, offset=0):
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
        init_pop = 'Initial population on reference',
        git_hash = 'VCS BASE repository version:'
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
    f.close()

    # Read table --- only read the first N columns, where N is the number of
    # column names found.
    data = pd.io.parsers.read_table(filename, sep='\s+', delim_whitespace=True,
               skiprows=start_line, skipfooter=skip_footer, names=column_names,
               comment='#')

    # Remove comment lines and convert all columns to numeric data.
    # Lines starting with a comment have been set to NaN in the iterations
    # column.
    try:
        data.dropna(subset=['iterations'], inplace=True)
    except TypeError:
        # Be slightly less efficient if using pandas version < 0.13.
        data = data.dropna(subset=['iterations'])
    data.reset_index(drop=True, inplace=True)
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
