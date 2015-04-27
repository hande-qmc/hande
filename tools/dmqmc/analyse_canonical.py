#!/usr/bin/env python

import pandas as pd
import os
import sys
script_dir = os.path.dirname(__file__)
sys.path.extend([os.path.join(script_dir, '../'),
                 os.path.join(script_dir, '../pyblock')])
import pyhande
import pyblock
import numpy as np


def main(filename):

    if len(filename) < 1:
        print("Usage: ./analyse_canonical.py files")
        sys.exit()

    (metadata, data) = pyhande.extract.extract_data_sets(filename)

    data.drop(labels='iterations', axis=1, inplace=True)
    data.rename(columns={'\\sum\\rho_HF_{ii}H_{ii}': 'num',
                '\\sum\\rho_HF_{ii}': 'denom'}, inplace=True)

    means = data.mean()
    covariances = data.cov()
    nsamples = len(data['num'])

    num = pd.DataFrame()
    num['mean'] = [means['num']]
    num['standard error'] = [np.sqrt(covariances['num']['num']/nsamples)]
    denom = pd.DataFrame()
    denom['mean'] = [means['denom']]
    denom['standard error'] = [np.sqrt(covariances['denom']['denom']/nsamples)]
    cov_thf = covariances['num']['denom']
    # The numerator and denominator are correlated for the
    # HF estimate for the total energy.
    e_thf = pyblock.error.ratio(num, denom, cov_thf, nsamples)
    e_thf.reset_index(inplace=True)

    results = pd.DataFrame()
    results['Beta'] = [b.split()[2].split(',')[0] for b in
                       metadata[0]['input'] if 'beta' in b]
    # E_0 and E_HF0 contain no denominator so the error is
    # just the standard error.
    results['E_0'] = [means['E_0']]
    results['E_0-Error'] = [np.sqrt(covariances['E_0']['E_0']/nsamples)]
    results['E_HF0'] = [means['E_HF0']]
    results['E_HF0-Error'] = [np.sqrt(covariances['E_HF0']['E_HF0']/nsamples)]
    results['E_THF'] = list(e_thf['mean'])
    results['E_THF-Error'] = list(e_thf['standard error'])

    print(results.to_string(index=False))


if __name__ == '__main__':

    main(sys.argv[1:])
