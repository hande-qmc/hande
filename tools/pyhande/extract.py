'''Data extraction from the (standard) output of a HANDE calculation.'''

import numpy
import pandas as pd

import json
import os
import re
import sys
import tempfile

import pyhande.legacy

def extract_data_sets(filenames):
    '''Extract QMC data tables from multiple HANDE calculations.

Parameters
----------
filenames : list of strings
    names of files containing HANDE QMC calculation output.

Returns
-------
data : list of (dict, :class:`pandas.DataFrame` or :class:`pandas.Series`)
    Calculation output represented by a tuple for each calculation, consisting
    of metadata (dict) and a :class:`pandas.DataFrame` (QMC calculations) or
    :class:`pandas.Series` (other calculations) containing the calculation
    output/results.

See Also
--------
``extract_data`` : underlying data extraction implementation.
'''

    data = []
    for filename in filenames:
        data.extend(extract_data(filename))
    return data

def extract_data(filename):
    '''Extract QMC data table from a HANDE calculation.

Parameters
----------
filename : string
    name of file containing the HANDE QMC calculation output.

Returns
-------
data_pairs : list of (dict, :class:`pandas.DataFrame` or :class:`pandas.Series`)
    Calculation output represented by a tuple for each calculation, consisting
    of metadata (dict) and a :class:`pandas.DataFrame` (QMC calculations) or
    :class:`pandas.Series` (other calculations) containing the calculation
    output/results.
'''

    # metadata from ...
    # ... output header
    md_generic_header = dict(
        UUID = 'Calculation UUID:',
        git_hash = 'VCS BASE repository version:',
    )
    # ... footer of output (ie after QMC data table)
    comms_footer = dict(
        min_psips_per_mpi_process = 'Min # of particles on a processor:',
        max_psips_per_mpi_process = 'Max # of particles on a processor:',
        mean_psips_per_mpi_process = 'Mean # of particles on a processor:',
        min_dets_per_mpi_process = 'Min # of determinants on a processor:',
        max_dets_per_mpi_process = 'Max # of determinants on a processor:',
        mean_dets_per_mpi_process = 'Mean # of determinants on a processor:',
        min_communication_time = 'Min time taken by walker communication:',
        max_communication_time = 'Max time taken by walker communication:',
        mean_communication_time = 'Mean time taken by walker communication:',
    )
    md_generic_footer = dict(
        wall_time = 'Wall time (seconds):',
        cpu_time = 'CPU time (per processor, seconds):'
    )

    data_pairs = []
    md_generic = {'input':[]}

    # RNG is for old calculations without proper calculation block headers.
    calc_block = re.compile('^ (FCI|FCIQMC|CCMC|DMQMC|Simple FCIQMC|'
                            'Hilbert space|Canonical energy|RNG)$')
    fci_block = re.compile('Exact|Lanczos|LAPACK|LANCZOS|RDM')

    # input block delimiters
    input_pattern = 'Input options'
    underline_regex = re.compile('----+')
    have_input = False

    calc_type = ''
    f = open(filename)
    for line in f:

        # Start of calculation block?
        match = calc_block.search(line)
        if match:
            calc_type = match.group().strip()

            if calc_type == 'FCI':
                metadata = _extract_json(f, find_start=True, max_end='subspace')
                metadata['calc_type'] = calc_type
            elif calc_type == 'Hilbert space':
                (metadata, data) = _extract_hilbert_data(f)
                metadata['calc_type'] = calc_type
                data_pairs.append((metadata, data))
            else:
                (metadata, data) = _extract_mc_calc(f)
                metadata['calc_type'] = calc_type
                data_pairs.append((metadata, data))
        elif calc_type == 'FCI' and fci_block.search(line):
            data = _extract_fci_data(f, line)
            data_pairs.append((metadata, data))
        elif not calc_type:
            # Generic header
            for (key, val) in md_generic_header.items():
                if val in line:
                    if key == 'git_hash':
                        md_generic[key] = next(f).strip()
                    else:
                        md_generic[key] = line.split()[-1]
            # Parse input block.
            if input_pattern in line:
                # skip next line and then start getting the input block.
                next(f)
                have_input = True
            elif have_input:
                if underline_regex.search(line):
                    have_input = False
                else:
                    md_generic['input'].append(line.strip())
        else:
            # QMC footer
            if calc_type != 'Hibert space' and calc_type != 'FCI':
                for (key, val) in comms_footer.items():
                    if val in line:
                        md_val = line.split()[-1].replace('s','')
                        data_pairs[-1][0][key] = float(md_val)
            # Generic footer
            for (key, val) in md_generic_footer.items():
                if val in line:
                    md_generic[key] = float(line.split()[-1])

    f.close()

    if data_pairs and 'system' not in data_pairs[0][0]:
        # Uhoh!  Have an old output with no JSON.  :-(
        # Note legacy metadata is *not* in the same format...
        md_legacy = pyhande.legacy.extract_metadata(filename)
        for (md, dat) in data_pairs:
            md.update(md_legacy)

    for (md, dat) in data_pairs:
        md.update(md_generic)

    return data_pairs

def _extract_mc_calc(fhandle):
    '''Extract metadata and calculation data for a QMC calculation.

Parameters
----------
fhandle : file
    File handle to a HANDE output file opened at the start of a QMC calculation
    block, as denoted by the section title.

Returns
-------
(metadata, data) : (dict, :class:`pandas.DataFrame`)
    Dictionary of calculation metadata (input values, defaults, etc) and the QMC
    data table obtained from the output file.
'''

    # Standard Monte Carlo table with an '# iterations ...' header.
    metadata = {}
    data = pd.DataFrame()
    for line in fhandle:
        if ' # iterations' in line:
            # Columns are separated by at least two spaces but each
            # column name can contain words separated by just one space.
            column_names = re.split('   *', line[3:].strip())
            # Work around pandas slow and very memory-hungry pure-python parser
            # by converting the data table into CSV format (and removing comments
            # whilst we're at it) and reading that in.
            tmp_csv = _convert_to_csv(fhandle)
            data = pd.io.parsers.read_csv(tmp_csv, names=column_names)
            os.remove(tmp_csv)
            # Done now -- return to main extraction procedure.
            break
        elif 'Start JSON block' in line:
            metadata = _extract_json(fhandle)

    if not data.empty:
        # If the number of iterations counter goes over 8 digits then the hande
        # output file prints stars.  This has now been fixed, however for legacy
        # reasons:
        for (i,iteration) in enumerate(data['iterations']):
            if numpy.isnan(iteration):
               data['iterations'][i] = \
                   i*metadata['qmc']['ncycles'] + data['iterations'][0]

        data = data.convert_objects(convert_numeric=True, copy=False)

        # Do we have an old table?  If so, rename the headings to the new
        # ones for convenience...
        data.rename(inplace=True, columns={
                        'Instant shift': 'Shift',
                        '\sum H_0j Nj': '\sum H_0j N_j',
                        '# D0': 'N_0',
                        '# particles': '# H psips',
                })

    return (metadata, data)

def _extract_fci_data(fhandle, title_line):
    '''Extract metadata and calculation data for an FCI calculation.

Parameters
----------
fhandle : file
    File handle to a HANDE output file opened at the start of a FCI calculation
    results block, as denoted by the section title.

.. note::

    The start should be at a block for the start of the results containing the
    eigenvalues.  A given FCI calculation can contain multiple such blocks.  The
    metadata is identical for all such blocks and is given at the start of the
    FCI section and hence should be extracted separately.

Returns
-------
data : :class:`pandas.Series`
    FCI eigenvalues obtained from the output file.
'''

    if 'Exact' in title_line or 'LAPACK' in title_line:
        title = 'FCI (LAPACK)'
    elif 'Lanczos' in title_line or 'LANCZOS' in title_line:
        title = 'FCI (Lanczos)'
    elif 'RDM' in title_line:
        title = 'FCI RDM'
    for line in fhandle:
        if 'State' in line:
            break
    eigvals = []
    for line in fhandle:
        if not line.split():
            break
        else:
            eigvals.append(float(line.split()[-1]))

    data = pd.Series(eigvals, index=list(range(1,len(eigvals)+1)))
    data.index.name = 'Eigenvalue'
    data.name = title

    return data

def _extract_hilbert_data(fhandle):
    '''Extract metadata and calculation data for a Hilbert space calculation.

Parameters
----------
fhandle : file
    File handle to a HANDE output file opened at the start of a Hilbert space
    calculation block, as denoted by the section title.

Returns
-------
(metadata, data) : (dict, :class:`pandas.Series`)
    Dictionary of calculation metadata (input values, defaults, etc) and the
    Hilbert space results obtained from the output file.
'''

    metadata = {}
    for line in fhandle:
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
        elif 'Start JSON block' in line:
            metadata = _extract_json(fhandle)
    return (metadata, data)

def _convert_to_csv(fhandle, comment='#'):
    '''Convert the HANDE data table in a file to CSV format.

Parameters
----------
fhandle: file
    python file handle of HANDE output file open at start of data table.
comment : string
    Single character which indicates a comment line if its the first
    non-whitespace character in a line.

.. note::

    We assume the table finishes at the next blank line.

Returns
-------
temp_file_name : string
    name of the CSV temporary file.
'''
    tmp = tempfile.NamedTemporaryFile(delete=False, mode='w')
    for line in fhandle:
        csv_line = ','.join(line.strip().split())
        if not csv_line:
            # blank line => end of data table.
            break
        elif not csv_line.startswith(comment):
            tmp.write(csv_line+'\n')
    tmp.close()
    return tmp.name

def _extract_json(fhandle, find_start=False, max_end=None):
    '''Extract JSON output from a HANDE output file.

Parameters
----------
fhandle : file
    File handle to a HANDE output file.
find_start : boolean
    If true, search for the start of the JSON block.  If false (default), then
    the file is assumed to be opened at the start of the JSON block.
max_end : string
    If find_start is True and max_end is not None, the search for the JSON block
    is aborted if a line containing max_end is found, in which case an empty
    dict is returned.

.. note::

    HANDE output contains blocks of output in JSON format.  The start of such
    a block is denoted by a line containing the string 'Start JSON block' and
    then end by a line containing the string 'End JSON block'.

Returns
-------
json_dict : dict
    JSON output loaded into a dictionary.
'''

    found_json = True
    if find_start:
        for line in fhandle:
            if 'Start JSON block' in line:
                break
            elif max_end is not None and max_end in line:
                found_json = False
                break
    json_text = ''
    if found_json:
        for line in fhandle:
            if ' End JSON block' in line:
                break
            else:
                json_text += line
    if json_text:
        return json.loads(json_text)
    else:
        return {}
