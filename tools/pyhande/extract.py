'''Data extraction from the (standard) output of a HANDE calculation.'''

import numpy
import pandas as pd
import os
import re
import sys
import tempfile

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

    data = []
    metadata = []
    for filename in filenames:
        # For now, ignore output from any other calculations.
        (calc_metadata, calc_data, junk) = extract_data(filename)
        data.append(calc_data)
        metadata.append(calc_metadata)
    return (pd.concat(metadata), pd.concat(data))

def extract_data(filename):
    '''Extract QMC data table from a HANDE calculation.

.. note::

    We assume that the data table is continuous (i.e. not split by anything
    other than comments) and that the final line of the table is within the
    final 2KB of the file.

Parameters
----------
filename : string
    name of file containing the HANDE QMC calculation output.

Returns
-------
metadata : :class:`pandas.Series`
    metadata (i.e. calculation information, parameters and settings) extracted
    from output file..
qmc_data : :class:`pandas.DataFrame`
    HANDE QMC data.
calc_data : list of `:class:`pandas.Series`
    data from other HANDE calculations (i.e. diagonalisation, Hilbert space
    estimation).
'''

    # metadata from ...

    # ... output header
    md_header = dict(
        UUID = 'Calculation UUID:',
        git_hash = 'VCS BASE repository version:',
    )

    # ... input (echoed in output)
    md_input = dict(
        calc_type = '^ *(simple_fciqmc|fciqmc|ccmc|ifciqmc|iccmc|dmqmc|idmqmc|\
                         |canonical_energy)',
        sym = r'\bsym\b +\d+',
        ms = r'\bms\b +-*\d+',
        nel = r'\bnel\b|\belectrons\b',
        tau = r'\btau\b',
        truncation = 'truncation_level',
        target = 'varyshift_target',
        shift_damping = 'shift_damping',
        mc_cycles = 'mc_cycles',
    )

    # ... main body of output (ie after input but before QMC data table)
    md_body = dict(
        nbasis = 'Number of basis functions:',
        seed = 'random number generator with a seed of',
        bit_length = 'Bit-length',
        ref = 'Reference determinant, |D0> = ',
        ref_energy = r'E0 = <D0|H|D0>',
        psingle = 'Probability of attempting a single',
        pdouble = 'Probability of attempting a double',
        init_pop = 'Initial population on',
    )

    # ... footer of output (ie after QMC data table)
    md_footer = dict(
        min_psips_per_mpi_process = 'Min # of particles on a processor:',
        max_psips_per_mpi_process = 'Max # of particles on a processor:',
        mean_psips_per_mpi_process = 'Mean # of particles on a processor:',
        min_dets_per_mpi_process = 'Min # of determinants on a processor:',
        max_dets_per_mpi_process = 'Max # of determinants on a processor:',
        mean_dets_per_mpi_process = 'Mean # of determinants on a processor:',
        min_communication_time = 'Min time taken by walker communication:',
        max_communication_time = 'Max time taken by walker communication:',
        mean_communication_time = 'Mean time taken by walker communication:',
        wall_time = 'Wall time (seconds):',
        cpu_time = 'CPU time (per processor, seconds):'
    )
    md_keys = [x for v in (md_header, md_input, md_body, md_footer)
            for x in v.keys()]
    md_int = ['sym', 'ms', 'nel', 'nbasis', 'truncation', 'seed', 'mc_cycles',
         'bit_length',  'min_dets_per_mpi_process', 'max_dets_per_mpi_process']

    md_float = ['tau', 'ref_energy', 'psingle', 'pdouble', 'init_pop',
         'min_psips_per_mpi_process', 'mean_psips_per_mpi_process',
         'mean_dets_per_mpi_process', 'mean_communication_time', 'wall_time',
         'cpu_time', 'min_communication_time', 'max_communication_time']

    input_pattern = 'Input options'
    underline_regex = re.compile('----+')

    # Read metadata and figure out how the start line of the data table.
    f = open(filename, 'r')
    start_line = 0
    metadata = pd.Series(index=md_keys, name='metadata')
    metadata['input'] = []
    unseen_calc = True
    have_git_hash_next = False
    have_input = 0
    calc_data = []
    calc_titles = ['diagonalisation results', 'RDM eigenvalues', 'Hilbert space']
    column_names = []
    for line in f:
        hit = False
        start_line += 1
        # extract metadata.
        if have_input <= 3:
            # input
            if have_input == 2 and line.strip() and \
                    not underline_regex.search(line):
                metadata['input'].append(line.strip())
                for (k,v) in md_input.items():
                    if re.search(v,line):
                        if k == 'calc_type':
                            if unseen_calc:
                                unseen_calc = False
                                metadata[k] = line.split()[0]
                        else:
                            metadata[k] = \
                                    extract_last_field(line, k, md_int, md_float)
                        hit = True
                        break
                if hit:
                    # only need to find each key once...
                    md_input.pop(k)
            if input_pattern in line or \
                    (have_input > 0 and underline_regex.search(line)):
                have_input += 1
        if have_input < 2:
            # header
            if have_git_hash_next:
                metadata['git_hash'] = line.split()[0]
                have_git_hash_next = False
            for (k,v) in md_header.items():
                if v in line:
                    if k == 'git_hash':
                        have_git_hash_next = True
                    elif k :
                        metadata[k] = extract_last_field(line, k, md_int, md_float)
        elif have_input > 3:
            # body
            hit = False
            for (k,v) in md_body.items():
                if v in line:
                    # Special cases for unusual formats...
                    if k == 'ref':
                        metadata[k] = ' '.join(line.split(v)).strip()
                    elif k == 'seed':
                        metadata[k] = int(float(line.split()[-1]))
                    else:
                        metadata[k] = extract_last_field(line, k, md_int, md_float)
                    hit = True
                    break
            if hit:
                md_body.pop(k)
            if any(k in line for k in calc_titles):
                (cdata, nread) = _extract_calc_data(f, line)
                calc_data.append(cdata)
                start_line += nread

        # Hunt for start of data table.
        if ' # iterations' in line:
            # Columns are separated by at least two spaces but each column name
            # can contain words separated by just one space.
            column_names = re.split('   *', line[3:].strip())
            break
    f.close()

    # Remapping to match old and new input files
    if metadata['calc_type'] == 'estimate_canonical_kinetic_energy':
        metadata['calc_type'] = 'kinetic_energy'

    if column_names:
        # Hunt for the end of the table.
        skip_footer = 0
        end_lines = _get_last_lines(filename)
        for line in end_lines[::-1]:
            if re.match('  *\d\d?', line):
                break
            else:
                skip_footer += 1

        # Extract meta data from the end of the calulation.
        for line in end_lines[-skip_footer:]:
            for (k,v) in md_footer.items():
                if v in line:
                    metadata[k] = extract_last_field(line, k, md_int, md_float)

        # Work around pandas slow and very memory-hungry pure-python parser
        # by converting the data table into CSV format (and stripping out
        # comments whilst we're at it) and reading that in.
        tmp_csv = _convert_to_csv(filename, start_line)
        qmc_data = pd.io.parsers.read_csv(tmp_csv, names=column_names)
        os.remove(tmp_csv)

        # If the number of iterations counter goes over 8 digits then the hande
        # output file prints stars.  This has now been fixed, however for
        # legacy reasons:
        for (i,iteration) in enumerate(qmc_data['iterations']):
            if numpy.isnan(iteration):
                qmc_data['iterations'][i] = \
                        i*metadata['mc_cycles'] + qmc_data['iterations'][0]

        qmc_data = qmc_data.convert_objects(convert_numeric=True, copy=False)

        # Do we have an old table?  If so, rename the headings to the new ones
        # for convenience...
        qmc_data.rename(inplace=True, columns={
                'Instant shift': 'Shift',
                '\sum H_0j Nj': '\sum H_0j N_j',
                '# D0': 'N_0',
                '# particles': '# H psips',
            })
    else:
        qmc_data = pd.DataFrame()

    return (metadata, qmc_data, calc_data)

def extract_last_field(line, key, md_int=None, md_float=None):
    '''Extract the final field from the last (space-separated) entry in a line.

Parameters
----------
line: string
    a line in the file
key: string
    the key to which has been matched
md_int: list
    a list of keys which should be stored as ints
md_float: list
    a list of keys which should be stored as floats

Returns
-------
val:
    The record, converted to an integer or a float if the key (the preceding
    fields in the line) is in md_int or md_float.
'''
    val = line.split()[-1]
    if val[-1] in ('.', ','):
        # Remove trailing commas and full-stops.
        val = val[:-1]
    if key in md_int:
        if "*" not in val:
            val = int(val)
    elif key in md_float:
        if val[-1] == 's':
            # Remove trailing s
            val = float(val[:-1])
        else:
            if "*" not in val:
                val = float(val)

    return val


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

def _extract_calc_data(fh, line):
    '''Extract non-QMC calculation output from a HANDE output file.

Parameters
----------

fh : file handle
    Open file handle to a HANDE output file, positioned at the start of a
    section containing the results from a non-QMC calculation.
line : string
    Contents of the current line in fh.

Returns
-------
data : `:class:`pandas.Series`
    Calculation data found.
nread : int
    Number of lines read.
'''
    nread = 0
    if any(key in line for key in ('Exact', 'Lanczos', 'LAPACK', 'LANCZOS', 'RDM')):
        if 'Exact' in line or 'LAPACK' in line:
            title = 'FCI (LAPACK)'
        elif 'Lanczos' in line or 'LANCZOS' in line:
            title = 'FCI (Lanczos)'
        elif 'RDM' in line:
            title = 'FCI RDM'
        for line in fh:
            nread += 1
            if 'State' in line:
                break
        eigvals = []
        for line in fh:
            nread += 1
            if not line.split():
                break
            eigvals.append(float(line.split()[-1]))
        data = pd.Series(eigvals, index=list(range(1,len(eigvals)+1)))
        data.index.name = 'Eigenvalue'
        data.name = title
    elif 'Hilbert space' in line:
        for line in fh:
            nread += 1
            if 'Monte-Carlo estimate of size of space is' in line:
                if '+/-' in line:
                    estimate = float(line.split()[-3])
                    std_err = float(line.split()[-1])
                    data = {'Monte Carlo estimate': estimate,
                            'standard error': std_err}
                else:
                    estimate = float(line.split()[-1])
                    data = {'Monte Carlo estimate': estimate}
                data = pd.Series(data)
                data.name = 'Hilbert space'
                break
            elif 'Size of space is' in line:
                # Deterministic value rather than MC estimate.
                data = pd.Series({'size': float(line.split()[-1])})
                data.name = 'Hilbert space'
                break
    return (data, nread)


def extract_input(metadata, variable):
    '''Extract input information about specific variable.

Parameters
----------
metadata : :class:`pandas.DataFrame`
    metadata (i.e. calculation information, parameters and settings) extracted
    from output files.
variable : string
    input variable whose value is desired.

Returns
-------
value : string
    value associated with variable in input file.
'''

    value = [x.split()[2].split(',')[0] for x in metadata['input'] if variable in x]
    return (value)
