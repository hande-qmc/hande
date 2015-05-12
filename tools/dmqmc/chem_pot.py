#!/usr/bin/env python

import sys
import scipy as sc
import scipy.optimize
import math
import numpy as np
import argparse
import pandas as pd


class UEGSystem:
    '''
    UEG system class. Slight overkill but this code is copied from a bigger
    script.

    Everything is measured in Hartree atomic units.

    Parameters
    ----------
    args : :class:`ArgumentParser'
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
        containing single particle eigenvalues, sorted in increasing value of
        kinetic energy.
    deg_e : list of lists
        compressed list of spval containing unique eigenvalues and their
        degeneracy.
    de : float
        epsilon value for comparing floats.
    '''
    def __init__(self, args):

        # Seitz radius.
        self.rs = args.rs
        # Number of electrons.
        self.ne = args.ne
        # Kinetic energy cut-off.
        self.ecut = args.ecutoff
        # Spin polarisation.
        self.pol = args.pol
        # Box Length.
        self.L = self.rs*(4*self.ne*sc.pi/3.)**(1/3.)
        # k-space grid spacing.
        self.kfac = 2*sc.pi/self.L
        # Fermi energy (inifinite systems).
        self.ef = 0.5*(9.0*sc.pi*self.pol/(4.0*self.rs**3.0))**(2./3.)
        # Single particle eigenvalues and corresponding kvectors
        self.spval = self.sp_energies(self.kfac, self.ecut)
        # Compress single particle eigenvalues by degeneracy.
        self.deg_e = compress_spval(self.spval)
        # epsilon value for comparison of floats.
        self.root_de = 1e-14

    def sp_energies(self, kfac, ecut):
        '''
        Calculate the allowed single particle eigenvalues which can fit in the
        sphere in kspace determined by ecut.

        Params
        ------
        kfac : float
            kspace grid spacing.
        ecut : float
            energy cutoff.

        Returns
        -------
        spval : list.
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


class MolSystem:
    '''
    Molecular system class.

    Parameters
    ----------
    args : :class:`ArgumentParser'
        where -
    ne : int
        number of electrons.
    pol : int
        spin polarisation = 2 for fully polarised case.
    beta : float
        inverse temperature.
    filename : string
        filename containing single-particle eigenvalues.
'''

    def __init__(self, args):

        # Number of electrons.
        self.ne = args.ne
        # Spin polarisation.
        self.pol = args.pol
        # Single particle eigenvalues.
        self.spval = self.sp_energies(args.filename[0])
        # Compress single particle eigenvalues by degeneracy.
        self.deg_e = compress_spval(self.spval)
        # epsilon value for comparison of floats.
        self.root_de = 1e-14

    def sp_energies(self, filename):

        data = pd.DataFrame()
        data = pd.read_csv(filename, delim_whitespace=True, skiprows=7, header=None)
        data.rename(columns={0:'a', 1:'b', 2:'c', 3:'d', 4:'e'}, inplace=True)
        speig = np.array(data.loc[(data.b > 0) & (data.c == 0) & (data.d == 0) & (data.e == 0),'a'])

        return sorted(speig)


def compress_spval(spval):
    '''
    Compress the single particle eigenvalues so that we only consider
    unique values which vastly speeds up the k-space summations required.
    [todo] - Look at more clever optimisations (stars).

    Params
    ------
    spval : list
        list containing single particle eigenvalues.

    Returns
    -------
    spval : list of lists.
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
            deg_e.append([j, eval1])
            i += j
            j = 1
        it += 1

    deg_e.append([j, eval1])

    return deg_e


def nav_sum(mu, ne, spval, beta, pol):
    '''Calculate average number of particles in the system.

Parameters
----------
mu : float
    chemical potential.
ne : int
    number of electrons.
spval : list of lists
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

    N = sum(degen/(np.exp(-beta*(mu-eigv))+1) for (degen, eigv) in spval)

    return (2.0/pol)*N


def nav_diff(mu, ne, spval, beta, pol):
    '''Calculate the difference between current estimate for N and the desired
    number of particles.  Needed as a separate function for use with Scipy root
    finding.

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


def chem_pot_sum(system, beta, guess):
    '''Find the correct value of the chemical potential which results in a
    system of ne particles at temperature beta.

Parameters
----------
system : class
    system parameters.
beta : float
    inverse temperature being considered.
guess : float
    initial guess for chemical potential.

Returns
-------
chemical potential : float
    the chemical potential.
'''

    return (sc.optimize.fsolve(nav_diff, guess, args=(system.ne,
                               system.deg_e, beta, system.pol))[0])


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

    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('calc', help='calculation type, either ueg or read_in')
    args = parser.parse_known_args()[0]
    print args
    if args.calc == 'ueg':
        parser.add_argument('rs', type=int, help='Wigner-Seitz radius.')
        parser.add_argument('ecutoff', type=float, help='Plane wave cutoff '
                            'in units of 0.5*(2\pi/L)**2.')
        parser.add_argument('-t', '--use--fermi', action='store_true',
                            dest='fermi_temperature', default=False,
                            help='Interpret input beta as Beta = T_F/T')
    elif args.calc == 'read_in':
        parser.add_argument('-f', '--filename', nargs='+', help='integral_file')
    else:
        print ("Calculation type unkown. Allowed values = ueg or read_in")

    parser.add_argument('ne', type=int, help='Number of electrons.')
    parser.add_argument('pol', type=int, help='Polarisation, = 1 for '
                        'unpolarised system and 2 for polarised system.')
    parser.add_argument('beta', type=float, help='Inverse temperature.')
    args = parser.parse_args()

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
        print("beta must be greater than zero.")
        sys.exit(1)
    else:
        if (args.calc == 'ueg'):
            system = UEGSystem(args)
            if args.fermi_temperature:
                args.beta = args.beta / system.ef
            chem_pot = chem_pot_sum(system, args.beta, system.ef)
        else:
            system = MolSystem(args)
            chem_pot = chem_pot_sum(system, args.beta, system.spval[system.ne-1])

    print(chem_pot)

if __name__ == '__main__':

    main(sys.argv[1:])
