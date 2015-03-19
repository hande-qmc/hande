'''Functions for obtaining data subsets to be tested by testcode'''

import pyhande.extract

def extract_test_data(fname):
    '''Extract data from a HANDE file and select a desired subset to be tested.

Parameters
----------
fname : string
    filename containing HANDE calculation output.

Returns
-------
output : dict
    A dictionary consisting of key, value pairs of calculation names and data
    subsets (each in a :class:`panda.DataFrame`).  Currently we select the first
    5 eigenvalues from FCI calculations, 5 equally spaced entries from QMC
    calculations and all data in other cases.
'''

    (metadata, qmc_data, other_calcs) = pyhande.extract.extract_data(fname)

    output = {}
    if not qmc_data.empty:
        qmc_data.rename(columns=lambda x: x.replace(' ', '_'), inplace=True)
        # Compare every 1/4 of the calculation...
        indxs = [0] + [int((i*len(qmc_data))/4)-1 for i in range(1,5)]
        test_data = qmc_data.ix[indxs]
        output[metadata[(0, 'calc_type')]] = test_data
    for calc in other_calcs:
        if 'FCI' in calc.name:
            # Compare at most the first 5 eigenvalues (just to be quick/minimise
            # test output/handle different orders of states within Lanczos):
            select_data = calc[:min(len(calc),5)].to_frame().T
            select_data.rename(columns=lambda x: 'eigv_%s' % (x,), inplace=True)
        elif 'Hilbert space' in calc.name:
            select_data = calc.to_frame().T
            select_data.rename(columns=lambda x: x.replace(' ', '_'), inplace=True)
        output[calc.name] = select_data

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
    for (key, dat) in extract_test_data(fname).items():
        dat_dict = dat.to_dict('list')
        for dat_key, dat_entry in dat_dict.items():
            if dat_key in test_data:
                test_data[dat_key].extend(dat_entry)
            else:
                test_data[dat_key] = dat_entry

    return test_data
