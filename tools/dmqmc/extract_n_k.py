#!/usr/bin/env python
'''Extract the momentum distribution from a analysed DMQMC simulation.'''

import pandas as pd
import numpy as np
import sys

if (len(sys.argv) < 2):
    print ("Usage: extract_n_k.py file bval")
    sys.exit()

bval = float(sys.argv[2])

data = pd.read_csv(sys.argv[1], sep=r'\s+').groupby('Beta').get_group(bval)

mom = [c for c in data.columns.values if 'n_' in c and '_error' not in c]
mome = [c for c in data.columns.values if 'n_' in c and '_error' in c]

vals = [float(c.split('_')[1]) for c in mom]

n_k = (data[mom].transpose()).values
n_k_error = (data[mome].transpose()).values
n_k_error[np.isnan(n_k_error)] = 0
frame = pd.DataFrame({'Beta': bval, 'k': vals, 'n_k': n_k.ravel(), 'n_k_error':
                      n_k_error.ravel()})

print (frame.to_string(index=False))
