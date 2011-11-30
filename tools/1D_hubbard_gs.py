#!/usr/bin/env python2
'''Calculate the energy per site of the ground-state wavefunction of a half-filled 1D Hubbard model in the thermodynamic limit.

Usage: 1D_hubbard_gs.ps U

where U is the ratio U/t which defines the system.

See:
E.H. Lieb, F.Y. Wu, Physica A: Statistical Mechanics and its Applications, 321 (2003) 1-27.
E.H. Lieb, F.Y. Wu, Phys. Rev. Lett. 20 (1968) 1445-1448.'''

import scipy
import scipy.special
import scipy.integrate
import sys

def wu_lieb_1d(U):
    '''Return the energy per site of a half-filled 1D Hubbard model in the thermodynamic limit, as given by Wu and Lieb.'''

    integrand = lambda w, U: (-4*scipy.special.j0(w)*scipy.special.j1(w))/(w*(1+scipy.exp(w*U/2)))

    (energy, err) = scipy.integrate.quad(integrand, 0, scipy.inf, args=(U))

    return (energy, err)

if __name__ == '__main__':

    if len(sys.argv) != 2:
        print __doc__
        sys.exit(1)

    U = float(sys.argv[1])

    (energy, err) = wu_lieb_1d(U)

    print 'Energy of the U=%dt half-filled 1D Hubbard model in the thermodynamic limit: %.16ft.\nEstimated absolute error in numerical integration: %g.'  % (U, energy, err)
