'''Attempt to remove population control bias by reweighting averages.'''

import pandas as pd
import numpy as np
from math import exp

def reweight(data, mc_cycles, tstep, weight_itts, mean_shift,
weight_key='Shift', arith_mean=False):
    '''Reweight using population control to reduce population control bias.
See C. J. Umirigar et. al. J. Chem. Phys. 99, 2865 (1993) equations 14 ..., 20. 

Rewieght estimators linear in the number of psips by a factor:
    W(tau, N) = \\Pi^{N-1}_{m=0} e^{-\\delta A \\tau S(\\tau - m)}

Where A is the number of steps per shift update cycle
S(\\tau - m) is the Shift on itteration \\tau - m,
\delta \\tau is the time step.
m is the number of itterations to reweight over.

Parameters
----------
data : :class:`pandas.DataFrame`
    HANDE QMC data.
tstep : float
    The time step used in the weight factor.
mc_cycles : int
    The number of monte carlo cycles per update step.
weight_itts: integer
    The number of iterations to reweight over.
mean_shift: float
    The mean shift (prevents weights becoming to big)
weight_key: string
    Column to generate the reweighting data.
geom_mean: bool
    Reweight using the geometric mean
Returns
-------
data : :class:`pandas.DataFrame`
    HANDE QMC data. with weights appended
'''
    weights = []
    to_prod = np.exp(-tstep*mc_cycles*(data[weight_key].values-mean_shift))
    if len(data[weight_key]) > 0:
        weights.append(to_prod[0])
    for i in range(1, min(weight_itts, len(data[weight_key]))):
        weights.append(weights[i-1]*to_prod[i])
    for i in range(weight_itts, len(data[weight_key])):
        weights.append(weights[i-1]*to_prod[i] / to_prod[i-weight_itts])
    if arith_mean:
        for i in range(weight_itts, len(data[weight_key])):
           weights[i] = weights[i]*arith_series(tstep, mc_cycles,
           data[weight_key].values[i], data[weight_key].values[i-weight_itts])

    data['Weight'] = weights

    return data

def arith_series(tstep, mc_cycles, weight_now, weight_before):
    # Handle division by zero!!!
    if weight_now != weight_before:
        series = (1/float(mc_cycles)) * ((1 - exp(-tstep*mc_cycles*(weight_before\
                  - weight_now)))/(1 - exp(-tstep*(weight_before - weight_now))))
    else:
        series = 1.0

    return series
