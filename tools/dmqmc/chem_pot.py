#!/usr/bin/env python

import sys
import scipy as sc
import scipy.optimize
import math
import numpy as np
import argparse

class System:
    '''
    UEG system class. Slight overkill but this code is copied from a bigger script.

    Everything is measured in Hartree atomic units.

    Parameters
    ----------
    args : list of command line arguments
        where -
    rs : float
        seitz radius.
    ne : int
        number of electrons.
    ecutoff : float
        plane wave cutoff.
    pol : int
        spin polarisation = 2 for fully polarised case.
    beta : float
        inverse temperature.

    Returns
    -------
    L : float
        box length.
    kfac : float
        kspace grid spacing.
    spval : numpy array
        containing single particle eigenvalues, sorted in increasing value of kinetic energy.
    deg_e : numpy array
        compressed list of spval containing unique eigenvalues and their degeneracy.
    de : float
        epsilon value for comparing floats.
    '''
    def __init__(self, args):

        # Seitz radius.
        self.rs = float(args[0])
        # Number of electrons.
        self.ne = int(args[1])
        # Kinetic energy cut-off.
        self.ecut = float(args[2])
        # Spin polarisation.
        self.pol = int(args[3])
        # Box Length.
        self.L = self.rs*(4*self.ne*sc.pi/3.)**(1/3.)
        # k-space grid spacing.
        self.kfac = 2*sc.pi/self.L
        # Fermi energy (inifinite systems).
        self.ef = 0.5*(9.0*sc.pi*self.pol/(4.0*self.rs**3.0))**(2./3.)
        # Single particle eigenvalues and corresponding kvectors
        self.spval = self.sp_energies(self.kfac, self.ecut)
        # Compress single particle eigenvalues by degeneracy.
        self.deg_e = self.compress_spval(self.spval)
        # epsilon value for comparison of floats.
        self.root_de = 1e-14

    def sp_energies(self, kfac, ecut):
        '''
        Calculate the allowed single particle eigenvalues which can fit in the sphere in kspace determined by ecut.

        Params
        ------
        kfac: float
            kspace grid spacing.
        ecut: float
            energy cutoff.

        Returns
        -------
        spval: numpy array.
            list containing single particle eigenvalues.
        '''

        # Scaled Units to match with HANDE.
        # So ecut is measured in units of 1/kfac^2.
        nmax = int(math.ceil(np.sqrt((2*ecut))))

        spval = []
        vec = []
        kval = []

        for ni in range(-nmax, nmax+1):
            for nj in range(-nmax, nmax+1):
                for nk in range(-nmax, nmax+1):
                    spe = 0.5*(ni**2 + nj**2 + nk**2)
                    if (spe <= ecut):
                        # Reintroduce 2 \pi / L factor.
                        spval.append(kfac**2*spe)

        # Sort in terms of increasing energy.
        spval.sort()

        return spval

    def compress_spval(self, spval):
        '''
        Compress the single particle eigenvalues so that we only consider unique
        values which vastly speeds up the k-space summations required.
        [todo] - Look at more clever optimisations (stars).

        Params
        ------
        spval: list
            list containing single particle eigenvalues.

        Returns
        -------
        spval: numpy array.
            list containing single particle eigenvalues.
        '''

        # Work out the degeneracy of each eigenvalue.
        j = 1
        i = 0
        it = 0
        deg_e = []
        deg_k = []

        while it < len(spval)-1:
            eval1 = spval[i]
            eval2 = spval[i+j]
            if eval2 == eval1:
                j += 1
            else:
                deg_e.append([j,eval1])
                i += j
                j = 1
            it += 1

        deg_e.append([j,eval1])

        return deg_e

def nav_sum(mu, ne, spval, beta, pol):
    '''Calculate average number of particles in the system.

Parameters
----------
mu : float
    chemical potential.
ne : int
    number of electrons.
spval : list of floats
    single particle energies and their degeneracies.
beta : float
    inverse temperature.
pol : int
    spin polarisation.

Returns
-------
N : float
    average number of particles.
'''

    N = 0

    for speig in range(0,len(spval)):

        N = N + spval[speig][0]/(np.exp(-beta*(mu-spval[speig][1]))+1)

    return (2.0/pol)*N

def nav_diff(mu, ne, spval, beta, pol):
    '''Calculate the difference between current estimate for N and the desired number of particles.
Needed as a separate function for use with Scipy root finding.

Parameters
----------
mu : float
    chemical potential.
ne : int
    number of electrons.
spval : list of floats
    single particle energies and their degeneracies.
beta : float
    inverse temperature.
pol : int
    spin polarisation.

Returns
-------
diff : float
    differenc between N and ne.
'''

    return nav_sum(mu, ne, spval, beta, pol) - ne

def chem_pot_sum(system, beta):
    '''Find the correct value of the chemical potential which results in a system of ne particles at temperature beta.

Parameters
----------
system : class
    system parameters.
beta : float
    inverse temperature being considered.

Returns
-------
chemical potential: float
    the chemical potential.
'''

    return sc.optimize.fsolve(nav_diff, system.ef, args=(system.ne, system.deg_e, beta, system.pol))[0]

def parse_args(args):
    '''Parse command-line arguments.

Parameters
----------
args : list of strings
    command-line arguments.

Returns
-------
args : :class:`ArgumentParser`
    Arguments read in from command line.
'''

    parser = argparse.ArgumentParser(usage = __doc__)
    parser.add_argument('rs', type=float, help='Wigner-Seitz radius.')
    parser.add_argument('ne', type=int, help='Number of electrons.')
    parser.add_argument('ecutoff', type=float, help='Plane wave cutoff in units of 0.5*(2\pi/L)**2.')
    parser.add_argument('pol', type=int, help='Polarisation, = 1 for unpolarised system and 2 for polarised system.')
    parser.add_argument('beta', type=float, help='Inverse temperature.')
    parser.add_argument('-t','--use--fermi', action='store_true', dest='fermi_temperature',
                      default=False, help='Interpret input beta as Beta = T_F/T')
    args = parser.parse_args(args)

    return args

def main(args):
    '''Calculate the chemical potential for the UEG.

Parameters
----------
args : list of numers.
    command-line arguments.

Returns
-------
None.
'''
    args = parse_args(args)
    if args.beta == 0:
        print "beta must be greater than zero."
        sys.exit(1)
    else:
        sys_pars = [args.rs, args.ne, args.ecutoff, args.pol]
        system = System(sys_pars)
        if args.fermi_temperature: args.beta = args.beta / system.ef
        chem_pot = chem_pot_sum(system, args.beta)

    print chem_pot

if __name__ == '__main__':

    main(sys.argv[1:])
