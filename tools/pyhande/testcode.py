'''Functions for obtaining data subsets to be tested by testcode'''

import pyhande.extract

def extract_test_data(fname, underscore=True):
    '''Extract data from a HANDE file and select a desired subset to be tested.

Parameters
----------
fname : string
    filename containing HANDE calculation output.
underscore : boolean
    if true, replace spaces in the key names with underscores.

Returns
-------
output : dict
    A dictionary consisting of key, value pairs of calculation names and data
    subsets (each in a :class:`panda.DataFrame`).  Currently we select the first
    5 eigenvalues from FCI calculations, 5 equally spaced entries from QMC
    calculations and all data in other cases.
'''

    hande_out = pyhande.extract.extract_data(fname)

    output = {}
    for (calc_ind, (metadata, data)) in enumerate(hande_out):
        # Label each calculation by the index of the calculation (in order
        # extracted from the test file, which is fixed) to handle the case where
        # multiple calculations of the same type are present in the same output
        # file.
        calc_ind = ' [%i]' % (calc_ind)
        if metadata['calc_type'] == 'FCI':
            # Compare at most the first 5 eigenvalues (just to be quick/minimise
            # test output/handle different orders of states within Lanczos):
            select_data = data[:min(len(data),5)].to_frame().T
            if underscore:
                select_data.rename(columns=lambda x: 'eigv_%s' % (x,), inplace=True)
            else:
                select_data.rename(columns=lambda x: 'eigv %s' % (x,), inplace=True)
            output[data.name + calc_ind] = select_data
        else:
            if underscore:
                data.rename(columns=lambda x: x.replace(' ', '_'), inplace=True)
            if metadata['calc_type'] == 'Hilbert space':
                # Compare all bar the first entry (which contains a legitimate
                # NaN in the standard error, which is not handled well by
                # testcode).
                indxs = range(len(data))[1:]
                test_data = data.ix[indxs]
                output[metadata['calc_type'] + calc_ind] = test_data
            elif len(data) > 10:
                # Compare every 1/4 of (non-trivial) QMC calculation...
                indxs = [0] + [int((i*len(data))/4)-1 for i in range(1,5)]
                test_data = data.ix[indxs]
                output[metadata['calc_type'] + calc_ind] = test_data
            else:
                output[metadata['calc_type'] + calc_ind] = data

    return output

def testcode_data(fname):
    '''Extract test data in the format required by testcode2.

Parameters
----------
fname : string
    filename containing HANDE calculation output.

Returns
-------
output : dict
    A dictionary consisting of key, value pairs of a data name and list of
    associated values.   This is essentially the output from
    ``extract_test_data`` combined into to a single dictionary.

See Also
--------

``extract_test_data`` : used to select the subset of data to be tested.
'''

    test_data = {}
    for (key, dat) in extract_test_data(fname, False).items():
        dat_dict = dat.to_dict('list')
        for dat_key, dat_entry in dat_dict.items():
            if dat_key in test_data:
                test_data[dat_key].extend(dat_entry)
            else:
                test_data[dat_key] = dat_entry

    return test_data
