'''Obtain metadata for legacy output files (i.e. not containing JSON blocks).'''
import re

def extract_metadata(filename):
    '''Extract metadata from a legacy output file.

.. note::

    The metadata is in the original format and **not** in the format obtained
    from the new JSON output.

Parameters
----------
filename : string
    Name of file containing HANDE output.

Returns
-------
metadata : dict
    Metadata extracted from analysing the (plain text) output and input options.
'''

    # metadata from ...
    # ... input (echoed in output)
    md_input = dict(
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

    md_int = ['sym', 'ms', 'nel', 'nbasis', 'truncation', 'seed', 'mc_cycles',
         'bit_length', 'mc_cycles']
    md_float = ['tau', 'ref_energy', 'psingle', 'pdouble', 'init_pop',
         'shift_damping', 'target']

    metadata = {}

    input_pattern = 'Input options'
    underline_regex = re.compile('----+')

    fh = open(filename)
    for line in fh:
        if input_pattern in line:
            next(fh)
            break
    # Parse metadata from input block.
    for line in fh:
        if underline_regex.search(line):
            break
        for (key, pattern) in md_input.items():
            if re.search(pattern, line):
                val = line.split()[-1]
                if val[-1] == ',':
                    val = val[:-1]
                if key in md_int:
                    metadata[key] = int(float(val))
                elif key in md_float:
                    metadata[key] = float(val)
                else:
                    metadata[key] = val
    # Parse metadata from body.
    for line in fh:
        if '# iterations' in line:
            # Finished with body! (Only worry about stopping early if there's
            # potentially lots of output left...
            break
        for (key, pattern) in md_body.items():
            if pattern in line:
                if key in md_int:
                    metadata[key] = int(float(line.split()[-1]))
                elif key in md_float:
                    metadata[key] = float(line.split()[-1])
                elif key == 'ref':
                    metadata[key] = ' '.join(line.split(pattern)).strip()
                else:
                    metadata[key] = line.split()[-1]

    return metadata

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
