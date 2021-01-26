'''Attempt to remove the population control bias by reweighting estimates.'''

import numpy as np

def reweight(data, mc_cycles, tstep, weight_history, mean_shift,
             weight_key='Shift'):
    '''Reweight using population control to reduce population control bias.

Reweight estimators linear in the number of psips by the factor:

.. math::

    W(\\tau, N) = \\Pi^{N-1}_{m=0} e^{-A \\delta \\tau S(\\tau - m\\delta\\tau)}

where :math:`A` is the number of steps per shift update cycle,
:math:`\\delta\\tau` is the time step and :math:`S(\\tau - m\\delta\\tau)` is
the shift at time :math:`\\tau - m\\delta\\tau`, and :math:`m` is the number of
iterations to reweight over.

See [Umrigar93]_ Eqs. 14-20 for details and [Vigor15]_ for use in FCIQMC.

Parameters
----------
data : :class:`pandas.DataFrame`
    HANDE QMC data.
tstep : float
    The time step used in the weight factor.
mc_cycles : int
    The number of monte carlo cycles per update step.
weight_history: integer
    The number of iterations to reweight over.
mean_shift: float
    The mean shift.  Used to prevent weights becoming too big.
weight_key: string
    Column to generate the reweighting data.

Returns
-------
weight : List[float]
    List of weights.
    

References
----------
Umrigar93
    C.J. Umirigar et al., J. Chem. Phys. 99, 2865 (1993)
Vigor15
    W.A. Vigor, et al., J. Chem. Phys. 142, 104101 (2015).
'''
    if tstep <= 0.0:
        raise ValueError("QMC time step tstep should be > 0 and not passed in "
                         f"value '{tstep}'.")
    if mc_cycles <= 0:
        raise ValueError("#Monte Carlo cycles mc_cycles should be > 0 and not"
                         f"passed in value '{mc_cycles}'.")

    weights = []
    to_prod = np.exp(-tstep*mc_cycles*(data[weight_key].values-mean_shift))
    if len(data[weight_key]) > 0:
        weights.append(to_prod[0])
    for i in range(1, min(weight_history, len(data[weight_key]))):
        weights.append(weights[i-1]*to_prod[i])
    for i in range(weight_history, len(data[weight_key])):
        weights.append(weights[i-1]*to_prod[i] / to_prod[i-weight_history])


    return weights

def arith_series(tstep, mc_cycles, weight_now, weight_before):
    
    if tstep <= 0.0:
        raise ValueError("QMC time step tstep should be > 0 and not passed in "
                         f"value '{tstep}'.")
    if mc_cycles <= 0:
        raise ValueError("#Monte Carlo cycles mc_cycles should be > 0 and not"
                         f"passed in value '{mc_cycles}'.")
    
    # Handle division by zero!!!
    if weight_now != weight_before:
        tdweight = tstep*(weight_before - weight_now)
        ratio = ((1 - exp(-mc_cycles*tdweight))/(1 - exp(-tdweight)))
        series = (1/float(mc_cycles)) * ratio
    else:
        series = 1.0

    return series