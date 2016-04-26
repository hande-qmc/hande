'''Extract data from the output of a HANDE calculation.

.. note::

    All :mod:`pyhande` analysis procedures assume data is in the format
    produced by :func:`extract_data` and :func:`extract_data_sets`.
'''

import numpy
import pandas as pd

import json
import os
import re
import sys
import tempfile
import warnings

import bz2
import gzip
try:
    import lzma
except ImportError:
    pass

import pyhande.legacy

def extract_data_sets(filenames):
    '''Extract QMC data tables from multiple HANDE calculations.

Parameters
----------
filenames : list of strings
    names of files containing HANDE QMC calculation output.

    .. note::

        Files compressed with gzip, bzip2 or xz (python 3 only) are
        automatically decompressed.

Returns
-------
data : list of (dict, :class:`pandas.DataFrame` or :class:`pandas.Series`)
    Calculation output represented by a tuple for each calculation, consisting
    of metadata (dict) and a :class:`pandas.DataFrame` (MC calculations) or
    :class:`pandas.Series` (other calculations) containing the calculation
    output/results.

See Also
--------
:func:`extract_data` : underlying data extraction implementation.
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

    .. note::

        Files compressed with gzip, bzip2 or xz (python 3 only) are
        automatically decompressed.

Returns
-------
data_pairs : list of (dict, :class:`pandas.DataFrame` or :class:`pandas.Series`)
    Calculation output represented by a tuple for each calculation, consisting
    of metadata (dict) and a :class:`pandas.DataFrame` (MC calculations) or
    :class:`pandas.Series` (other calculations) containing the calculation
    output/results.
'''

    # metadata from ...
    # ... output header
    md_generic_header = dict(
        UUID = 'Calculation UUID:',
        git_hash = re.compile('git sha1 hash:|VCS BASE repository version:'),
        hande_version = 'HANDE version:',
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

    calc_types = '(FCI|FCIQMC|CCMC|DMQMC|Simple FCIQMC|'\
                            'Hilbert space|Canonical energy)'
    calc_block = re.compile('^ '+calc_types+'$')
    fci_block = re.compile('Exact|Lanczos|LAPACK|LANCZOS|RDM')
    timing_pattern = re.compile('^ '+calc_types+
                                ' (?:estimation|calculation) *: ([0-9.]+)$')

    # input block delimiters
    input_pattern = 'Input options'
    underline_regex = re.compile('----+')
    have_input = False

    calc_type = ''
    timings = []
    (f, compressed) = _open_file(filename)
    for line in f:

        # Start of calculation block?
        match = calc_block.search(line)
        if 'RNG' in line:
            warnings.warn('Data extraction identified from an RNG block no '
	                  'longer supported.  Please use an older version of '
			  'pyhande.')
        elif match:
            calc_type = match.group().strip()

            if calc_type == 'FCI':
                metadata = _extract_json(f, find_start=True, max_end='subspace')
                metadata['calc_type'] = calc_type
            elif calc_type == 'Hilbert space':
                (metadata, data) = _extract_hilbert_data(f)
                metadata['calc_type'] = calc_type
                data_pairs.append((metadata, data))
            else:
                (metadata, data, comment_data) = _extract_mc_calc(f, calc_type)
                metadata['calc_type'] = calc_type
                data_pairs.append((metadata, data))
                if calc_type == 'DMQMC' and not comment_data.empty:
                    # Also got some results in the comment_file...
                    metadata_rdm = metadata.copy()
                    metadata_rdm['calc_type'] = 'DMQMC (RDM)'
                    data_pairs.append((metadata_rdm, comment_data))
        elif calc_type == 'FCI' and fci_block.search(line):
            data = _extract_fci_data(f, line)
            data_pairs.append((metadata, data))
        elif not calc_type:
            # Generic header
            for (key, val) in md_generic_header.items():
                if key == 'git_hash':
                    if val.search(line):
                        md_generic[key] = next(f).strip()
                elif val in line:
                    md_generic[key] = line.split()[-1].strip('.')
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
                        # Check if Fortran has starred out the number
                        if "*" in md_val:
                            data_pairs[-1][0][key] = float('nan')
                        else:
                            data_pairs[-1][0][key] = float(md_val)
            # Timing summary
            if 'Timing breakdown' in line:
                # Skip underline and blank line
                next(f)
                next(f)
                next(f)
                next(f)
                for line in f:
                    if not line.strip():
                        break
                    m = timing_pattern.match(line)
                    if m: 
                        timings.append((m.group(1), float(m.group(2))))
            # Generic footer
            for (key, val) in md_generic_footer.items():
                if val in line:
                    md_generic[key] = float(line.split()[-1])

    f.close()

    if data_pairs and 'system' not in data_pairs[0][0]:
        # Uhoh!  Have an old output with no JSON.  :-(
        # Note legacy metadata is *not* in the same format...
        (f, compressed) = _open_file(filename)
        md_legacy = pyhande.legacy.extract_metadata(f)
        for (md, dat) in data_pairs:
            md.update(md_legacy)
        f.close()

    for (md, dat) in data_pairs:
        md.update(md_generic)
    if timings:
        # Assume order calculations are run in is the same as the timing report.
        for ((md, dat), (calc_type, time)) in zip(data_pairs, timings):
            md.update({'calculation_time':time})

    return data_pairs

def _extract_mc_calc(fhandle, calc_type):
    '''Extract metadata and calculation data for a QMC calculation.

Parameters
----------
fhandle : file
    File handle to a HANDE output file opened at the start of a QMC calculation
    block, as denoted by the section title.
calc_type : string
    Type of calculation being analysed, e.g. 'CCMC', 'DMQMC' or 'FCIQMC'.
    Currently only used in dealing with data encoded in comment lines.

Returns
-------
(metadata, data, comment_data) : (dict, :class:`pandas.DataFrame`, :class:`pandas.DataFrame`)
    Dictionary of calculation metadata (input values, defaults, etc), the QMC
    data table obtained from the output file and any data extracted from the
    comment lines.
'''

    # Standard Monte Carlo table with an '# iterations ...' header.
    metadata = {}
    data = pd.DataFrame()
    comment_data = None
    data_csv = None
    for line in fhandle:
        if ' # iterations' in line:
            # Columns are separated by at least two spaces but each
            # column name can contain words separated by just one space.
            column_names = re.split('   *', line[3:].strip())
            # Work around pandas slow and very memory-hungry pure-python parser
            # by converting the data table into CSV format (and removing
            # comment_file whilst we're at it) and reading that in.
            (data_csv, comment_file) = _convert_to_csv(fhandle)
            # Done now -- return to main extraction procedure.
            break
        elif 'Start JSON block' in line:
            metadata = _extract_json(fhandle)
    # It's possible that the output file didn't have any info, so we need to
    # test data_csv
    if data_csv:
        data = pd.io.parsers.read_csv(data_csv, names=column_names)
        if calc_type == 'DMQMC':
            comment_data = _extract_dmqmc_data(comment_file)

        os.remove(data_csv)
        os.remove(comment_file)

    if not data.empty:
        # If the number of iterations counter goes over 8 digits then the hande
        # output file prints stars.  This has now been fixed, however for legacy
        # reasons:
        data['iterations'].replace('\*+', -1, regex=True, inplace=True)
        try:
            data['iterations'] = pd.to_numeric(data['iterations'])
        except AttributeError:
            data = data.convert_objects(convert_numeric=True, copy=False)
        for (i,iteration) in data['iterations'].iteritems():
            if iteration < 0:
                data.loc[i,'iterations'] = \
                   i*metadata['qmc']['ncycles'] + data.loc[0,'iterations']

        # Do we have an old table?  If so, rename the headings to the new
        # ones for convenience...
        data.rename(inplace=True, columns={
                        'Instant shift': 'Shift',
                        '\sum H_0j Nj': '\sum H_0j N_j',
                        '# D0': 'N_0',
                        '# particles': '# H psips',
                })

    return (metadata, data, comment_data)

def _extract_dmqmc_data(comment_file):
    '''Extract data from comments produced by a DMQMC calculation.

Parameters
----------
comment_file : string
    Name of file containing the comment lines, as produced by
    :func:`_convert_to_csv`.

Returns
-------
data : :class:`pandas.DataFrame`
    Data contained in the comments (e.g. ground state RDM quantities).
'''
    data = {}
    keys = ['RDM trace', 'von Neumann', 'concurrence']
    f = open(comment_file)
    for line in f:
        for key in keys:
            if key in line:
                val = float(line.split()[-1])
                if key in data:
                    data[key].append(val)
                else:
                    data[key] = [val]
    f.close()
    data = pd.DataFrame(data)
    data.index.name = 'beta loop'
    return data

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
(metadata, data) : (dict, :class:`pandas.DataFrame`)
    Dictionary of calculation metadata (input values, defaults, etc) and the
    Hilbert space results obtained from the output file.
'''

    metadata = {}
    for line in fhandle:
        if '# iterations' in line:
            # Columns are separated by at least two spaces but each
            # column name can contain words separated by just one space.
            column_names = re.split('   *', line[3:].strip())
            (data_csv, junk) = _convert_to_csv(fhandle, parse_comments=False)
            data = pd.io.parsers.read_csv(data_csv, names=column_names)
            os.remove(data_csv)
            data_table = True
            break
        elif 'Monte-Carlo estimate of size of space is' in line:
            # Old Monte Carlo algorithm (no iterations, std. err. estimate from
            # different MPI ranks).
            if '+/-' in line:
                estimate = float(line.split()[-3])
                std_err = float(line.split()[-1])
            else:
                estimate = float(line.split()[-1])
                std_err = float('nan')
            data_table = False
            break
        elif 'Size of space is' in line:
            # Deterministic value rather than MC estimate.
            estimate = float(line.split()[-1])
            std_err = 0.0
            data_table = False
            break
        elif 'Start JSON block' in line:
            metadata = _extract_json(fhandle)
    if not data_table:
        data = {'iterations': [1],
                'space size': [estimate],
                'mean': [estimate],
                'std. err.': [std_err]}
        data = pd.DataFrame(data)
    data.name = 'Hilbert space'
    return (metadata, data)

def _convert_to_csv(fhandle, comment='#', parse_comments=True):
    '''Convert the HANDE data table in a file to CSV format.

Parameters
----------
fhandle: file
    python file handle of HANDE output file open at start of data table.
comment : string
    Single character which indicates a comment line if its the first
    non-whitespace character in a line.
parse_comments : boolean
    If true , also save the comment lines to a separate file.

.. note::

    We assume the table finishes at the next blank line.

Returns
-------
temp_filename : string
    name of the CSV temporary file containing the data table.
comment_file_filename : string
    name of file containing the comment lines extracted from the data table
    and an empty string if parse_comments.if False.
'''
    data = tempfile.NamedTemporaryFile(delete=False, mode='w')
    if parse_comments:
        comment_file = tempfile.NamedTemporaryFile(delete=False, mode='w')
        comment_fname = comment_file.name
    else:
        comment_fname = ''
    for line in fhandle:
        csv_line = ','.join(line.strip().split())
        if not csv_line:
            # blank line => end of data table.
            break
        elif csv_line.startswith(comment):
            if parse_comments:
                comment_file.write(line)
        else:
            data.write(csv_line+'\n')
    data.close()
    if parse_comments:
        comment_file.close()
    return (data.name, comment_fname)

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

def _open_file(fname):
    '''Open a gzip/bz2/ascii file.

Parameters
----------
fname : string
    name of file to be opened.

Returns
-------
(fhandle, compressed) : (file, bool)
    fhandle is the handle of opened file (a BZ2File or GzipFile or file object,
    depending upon the file passed in) and compressed is set to True if
    a compressed file (of any type) is detected and False otherwise.
'''

    # See http://stackoverflow.com/a/13044946/3412233.
    py3 = sys.version_info[0] == 3
    if py3:
        magic_nos = {
            b'\x1f\x8b\x08': gzip.open,             # gzip
            b'\x42\x5a\x68': bz2.open,              # bz2
            b'\xfd\x37\x7a\x58\x5a\x00': lzma.open, # xz
        }
    else:
        magic_nos = {
            b'\x1f\x8b\x08': gzip.GzipFile,      # gzip
            b'\x42\x5a\x68': bz2.BZ2File,        # bz2
            b'\xfd\x37\x7a\x58\x5a\x00': 'lzma', # xz
        }

    max_magic_len = max(len(x) for x in magic_nos)
    f = open(fname, 'rb')
    file_start = f.read(max_magic_len)
    f.close()
    compressed = False
    for (magic, fzopen) in magic_nos.items():
        if file_start.startswith(magic):
            compressed = True
            if py3:
                f = fzopen(fname, encoding='utf-8', mode='rt')
            elif fzopen == 'lzma':
                raise RuntimeError('xz compression not available on python 2.')
            else:
                f = fzopen(fname)
    if not compressed:
        f = open(fname)
    return (f, compressed)
