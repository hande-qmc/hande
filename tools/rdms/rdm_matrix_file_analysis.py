#!/usr/bin/env python
''' rdm_energy_estimate.py [options] RDM_1 RDM_2 ... RDM_N

Average together RDM files produced during an FCIQMC calculation and
perform additional analysis such as calculating the energy estimate from
the RDM profile.
'''

import os
import sys
import pkgutil
import argparse
import pandas as pd

try:
    from pyhande.extract import extract_data
    from pyblock.error import ratio
except ModuleNotFoundError:
    _script_dir = os.path.dirname(os.path.abspath(__file__))
    if not pkgutil.find_loader('pyhande'):
        sys.path.append(os.path.join(_script_dir, '../pyhande'))
    from pyhande.extract import extract_data
    from pyblock.error import ratio


def parse_arguments(arguments):
    ''' Parse command-line arguments.

    Parameters
    ----------
    arguments : list of strings
        command-line arguments.

    Returns
    -------
    list_of_filenames : list of strings
        list HANDE RDM files from an FCIQMC calculation.
    options : :class:`ArgumentParser`
        Options read in from command line.

    Raises
    ------
    UserWarning
    '''

    parser = argparse.ArgumentParser(usage=__doc__)

    parser.add_argument('-mf', '--metadata_file', action='store', default=None,
                        type=str, dest='metadata_file', help='Define a file '
                        'name for which metadata related to the calculation '
                        'is stored. The file can be used to set the CAS,  '
                        'reference determinant, int_file, and complex rather '
                        'than manually providing them. Note, if the script is '
                        'run on a directory, we assume the int_file has the '
                        'same path as the RDM file(s).')
    parser.add_argument('-intf', '--int_file', action='store', default=None,
                        type=str, dest='int_file', help='Define a file '
                        'which contains the integrals for the system. ')
    parser.add_argument('-j', '--complex-integrals', action='store_true',
                        dest='complex', default=False, help='Specify the '
                        'symmetry of the orbitals. Only alters the symmetries '
                        'for finding integrals and RDM elements. Currently '
                        'there is no routines for addressing integrals with '
                        'an imaginary component.')
    parser.add_argument('-cas', '--cas', nargs=2, default=None, type=int,
                        help='Give the active space used during the '
                        'calculation, as space separated values: M N. Note '
                        'this is required for accurate analysis if a CAS was '
                        'used during the calculation.')
    parser.add_argument('-det', '--det', nargs='+', type=int, help='Provide '
                        'the occupied orbitals used to define the reference '
                        'of the system. Required when using an occupation '
                        'vector different than the auto-filled N orbitals.')
    parser.add_argument('filenames', nargs='+', help='HANDE files to analyse.')
    parser.parse_args(args=None if arguments else ['--help'])

    options = parser.parse_args(arguments)

    return options.filenames, options


def read_dump_file(dump):
    ''' A very basic function for reading and storing the integrals from an
    FCIDUMP. We assume that the format for the FCIDUMP file follows that of:
    https://hande.readthedocs.io/en/latest/manual/integrals.html#fcidump-format

    Warning we just store each spin symmetry so this can be very memory hungry.

    Parameters
    ----------
    dump : str
        The filename of the FCIDUMP to read in
        system information and integrals.

    Returns
    -------
    integrals : dictionary
        The integrals for the system, we store each spin symmetry in
        physicists notation.
    nel : int
        The number of electrons in the system.
    norb : int
        The number of spin orbitals in the system.
    uhf : boolean
        True if the integrals correspond to an unrestricted Hartree--Fock
        calculations.

    Raises
    ------
    RuntimeError
        If the integrals are found to be in the form z = x + iy.
    RuntimeError
        If the FCIDUMP file has the integrals and orbital indices contained
        on multiple lines.
    RuntimeError
        If one of nel or norb where not identified.
    '''

    integrals = {}
    nel = None
    norb = None
    uhf = False
    header_complete = False

    with open(dump, 'r') as stream:

        for line in stream:
            line = line.upper().replace('NELEC', 'NEL')

            if 'END' in line:
                header_complete = True
                if not uhf:
                    norb += norb
                    rhf_fac = 2
                else:
                    rhf_fac = 1
            elif not header_complete:
                line = 'LINE_PAD ' + line.replace('=', '')
                if 'NEL' in line:
                    nel = int(line.split('NEL')[1].split(',')[0])
                if 'NORB' in line:
                    norb = int(line.split('NORB')[1].split(',')[0])
                if 'UHF' in line:
                    uhf = 'TRUE' in line
            elif header_complete:
                if '(' in line and ')' in line:
                    raise RuntimeError('Complex integral values are not '
                                       'supported, please send patches!')

                line_data = line.split()

                if len(line_data) != 5:
                    raise RuntimeError('Format of FCIDUMP file is not '
                                       'recognized, please read the '
                                       'documentation for a format outline '
                                       'or send patches!')

                i, a, j, b = (int(rhf_fac*int(k)) for k in line_data[1:])
                integrals[i, j, a, b] = float(line_data[0])

                if uhf:
                    continue
                elif i != 0 and a != 0 and j != 0 and b != 0:
                    integrals[i-1, j-1, a-1, b-1] = integrals[i, j, a, b]
                    integrals[i, j-1, a, b-1] = integrals[i, j, a, b]
                    integrals[i-1, j, a-1, b] = integrals[i, j, a, b]
                elif i != 0 and a != 0:
                    integrals[i-1, j, a-1, b] = integrals[i, j, a, b]
                elif i != 0:
                    integrals[i-1, j, a, b] = integrals[i, j, a, b]

    if nel is None or norb is None:
        raise RuntimeError('The number of electrons or number of orbitals '
                           'was not found in the FCIDUMP, send patches!')

    return integrals, nel, norb, uhf


def read_density_matrix_file(rdm):
    ''' Read in the RDM file produces from an FCIQMC calculation and store
    the values in a dictionary indexed by orbitals.

    Parameters
    ----------
    rdm : str
        A file containing the RDM values.

    Returns
    -------
    gamma : dictionary
        A dictionary with the RDM values with the spin orbitals as the keys.
    '''
    gamma = {}
    with open(rdm, 'r') as stream:
        for line in stream:
            line_data = line.split()
            i, j, a, b = (int(k) for k in line_data[:-1])
            gamma_ijab = float(line_data[-1])
            gamma[i, j, a, b] = gamma_ijab

    return gamma


def orbital_index_from_pair_index(index_pair):
    ''' A python version of the orbs_from_index function within HANDE.
    For more information, see `orbs_from_index` in: lib/local/utils.F90
    We always return the smallest index first.

    Parameters
    ----------
    index_pair : int
        A single index representing a pair of orbitals.

    Returns
    -------
    index_one : int
        The first index of the pair of orbitals.
    index_two : int
        The second index of the pair of orbitals.
    '''
    index_one = int(1.5 + (2*index_pair-1.75)**0.5)
    index_two = index_pair - int(int((index_one-1)*(index_one-2))/2)

    if index_one > index_two:
        index_one, index_two = index_two, index_one

    return index_one, index_two


def get_two_particle(integrals, i, j, a, b, bcomplex):
    ''' Permute through the orbital symmetries for the given system to
    identify a non-zero value within the given array. If no value is found
    to be within, zero is returned.

    For the permutational symmetries see:
    1) http://vergil.chemistry.gatech.edu/notes/permsymm/permsymm.html
    2) `read_in_integrals` in: src/read_in.F90

    Parameters
    ----------
    integrals : dictionary
        The tensor array being searched for non-zero values given some
        pair of spin orbitals.
    i, j, a, b : int
        The spin orbitals
    bcomplex : boolean
        True if the integrals have a 4-permutation symmetry.

    Returns
    -------
    float
        Zero if no value is found, otherwise the non-zero value
        corresponding to the 4 orbitals provided is returned.
    '''
    if not bcomplex:
        symmetries = [
                [i, j, a, b],
                [j, i, b, a],
                [a, b, i, j],
                [b, a, j, i],
                [a, j, i, b],
                [b, i, j, a],
                [i, b, a, j],
                [j, a, b, i],
            ]
    else:
        symmetries = [
                [i, j, a, b],
                [j, i, b, a],
                [a, b, i, j],
                [b, a, j, i],
            ]

    for ii, jj, aa, bb in symmetries:
        if (ii, jj, aa, bb) in integrals.keys():
            return integrals[ii, jj, aa, bb]

    return 0.0


def search_twofold(integrals, i, a):
    ''' Search the rank-4 dictionary for two free orbital indices corresponding
    to a non-zero value.

    Parameters
    ----------
    integrals : dictionary
        A rank 4 dictionary containing the integral values for the system.
    i, a : int
        Orbital indices to permute.

    Returns
    -------
    float
        Zero if no value is found, otherwise the non-zero value
        corresponding to the 2 orbitals provided is returned.
    '''

    if (i, 0, a, 0) in integrals.keys():
        return integrals[i, 0, a, 0]
    elif (a, 0, i, 0) in integrals.keys():
        return integrals[a, 0, i, 0]

    return 0.0


def get_one_particle(integrals, a, b, Nfrozen, bcomplex):
    ''' Finds a non-zero one-particle integral within a rank 4 dictionary.
    If a CAS space is being used, calculate the appropriate one-particle
    integral given the CAS space.

    For more information on the one-particle calculation given a CAS space,
    see comments within `read_in_integrals` in: src/read_in.F90

    Parameters
    ----------
    integrals : dictionary
        A rank 4 dictionary containing the integral values for the system.
    a, b : int
        Orbital indices corresponding to the one-particle integral.
    Nfrozen : int
        The number of frozen orbitals for the CAS, zero if not CAS is being
        used in the calculation.
    bcomplex : boolean
        True if the integrals have a 4-permutation symmetry.

    Returns
    -------
    hab : float
        The one-particle Hamiltonian element.
    '''

    hab = search_twofold(integrals, a, b)

    if Nfrozen == 0:
        return hab

    for i in range(1, Nfrozen + 1):
        hab += get_two_particle(integrals, i, a, i, b, bcomplex)
        hab -= get_two_particle(integrals, i, a, b, i, bcomplex)

    return hab


def calculate_density_matrix_energy(gamma, integrals, nel, norb, Mfrozen,
                                    det, bcomplex):
    ''' Calculate the numerator and denominator of the energy expression
    from the RDM file. This is simply a python implementation of the
    `calc_rdm_energy` function in: src/replica_rdm.F90

    Parameters
    ----------
    gamma : dictionary
        A dictionary with the RDM values with the spin orbitals as the keys.
    integrals : dictionary
        A rank 4 dictionary containing the integral values for the system.
    nel : int
        The number of electrons in the system.
    norb : int
        The number of spin orbitals in the system.
    Nfrozen : int
        The number of frozen orbitals for the CAS, zero if not CAS is being
        used in the calculation.
    det : list of int
        The occupied spin orbitals for the system.
    bcomplex : boolean
        True if the integrals have a 4-permutation symmetry.

    Returns
    -------
    numerator : float
        The numerator for the energy expression from the RDM.
    trace : float
        The denominator for the energy expression from the RDM.
    '''

    numerator = 0.0
    trace = 0.0

    for ij in range(1, norb + 1):
        for ab in range(1, norb + 1):

            i, j = orbital_index_from_pair_index(ij)
            a, b = orbital_index_from_pair_index(ab)

            gamma_ij = get_two_particle(gamma, i, j, a, b, bcomplex)

            if gamma_ij == 0.0:
                continue
            elif ij == ab:
                trace += gamma_ij

            i, j, a, b = i+Mfrozen, j+Mfrozen, a+Mfrozen, b+Mfrozen

            v_ijab = get_two_particle(integrals, i, j, a, b, bcomplex)
            v_ijba = get_two_particle(integrals, i, j, b, a, bcomplex)
            vbar = v_ijab - v_ijba

            numerator += gamma_ij * vbar

            if i == a:
                h_jb = get_one_particle(integrals, j, b, Mfrozen, bcomplex)
                numerator += gamma_ij * h_jb/(nel-1)
            if j == b:
                h_ia = get_one_particle(integrals, i, a, Mfrozen, bcomplex)
                numerator += gamma_ij * h_ia/(nel-1)
            if i == b:
                h_ja = get_one_particle(integrals, j, a, Mfrozen, bcomplex)
                numerator -= gamma_ij * h_ja/(nel-1)
            if j == a:
                h_ib = get_one_particle(integrals, i, b, Mfrozen, bcomplex)
                numerator -= gamma_ij * h_ib/(nel-1)

    numerator = numerator*nel*(nel-1)*0.5

    core = get_two_particle(integrals, 0, 0, 0, 0, bcomplex)

    for i in range(1, Mfrozen + 1):
        core += get_one_particle(integrals, i, i, 0, bcomplex)

        for j in range(i + 1, Mfrozen + 1):
            v_ijab = get_two_particle(integrals, i, j, i, j, bcomplex)
            v_ijba = get_two_particle(integrals, i, j, j, i, bcomplex)
            vbar = v_ijab - v_ijba

            core += vbar

    reference = get_two_particle(integrals, 0, 0, 0, 0, bcomplex)

    occs = []
    for ifrzn in range(1, Mfrozen + 1):
        occs.append(ifrzn)
    for iactv in det:
        occs.append(Mfrozen + iactv)

    for indx, i in enumerate(occs):

        reference += get_one_particle(integrals, i, i, 0, bcomplex)

        for j in occs[indx+1:]:

            v_ijab = get_two_particle(integrals, i, j, i, j, bcomplex)
            v_ijba = get_two_particle(integrals, i, j, j, i, bcomplex)
            vbar = v_ijab - v_ijba

            reference += vbar

    numerator += (core - reference)*trace

    return numerator, trace


def stdout_average_gamma(gamma):
    ''' Calculate the average RDM matrix elements and the standard error
    given a data frame containing matrix elements sampled in different
    FCIQMC simulations. Report the average and error along with the spin
    orbital indices for the RDM matrix element to stdout.

    Parameters
    ----------
    gamma : :class:`pandas.DatFrame`
        A data frame containing the 4 spin orbitals labeling the
        matrix elements of the RDM, along with the matrix element sampled
        in an FCIQMC simulation.

    Returns
    -------
    None.
    '''
    grp = gamma.groupby('ijab')
    out = grp.mean()
    out['gamma_error'] = grp.sem()['gamma']

    print()
    print('Averaged RDM profile:')
    print(f' {"i":>5} {"j":>5} {"a":>5} {"b":>5} '
          f'{"gamma_kl":>24} {"gamma_kl_error":>24}')
    for ijab in out.index:
        ave, err = out.loc[ijab, 'gamma'], out.loc[ijab, 'gamma_error']
        i, j, a, b = ijab.split()
        print(f' {i:>5} {j:>5} {a:>5} {b:>5} {ave:> 24.12E} {err:>24.12E}')


def stdout_average_energy(estimates):
    ''' Calculate the average energy estimate from the RDM files produced
    during an FCIQMC simulation. For more information, see the DMQMC analysis
    function `analyse_observables` within: tools/pyhande/pyhande/dmqmc.py.
    Report the final energy by printing to stdout.

    Parameters
    ----------
    estimates : :class:`pandas.DatFrame`
        A data frame containing the numerators and denominators calculated
        with the read in RDM files.

    Returns
    -------
    None.
    '''
    df = pd.DataFrame(estimates)
    ave = df.mean()
    sem = df.sem()
    cov = df.cov()
    N = df.count()['trace']

    num = pd.DataFrame()
    num['mean'] = [ave['numerator']]
    num['standard error'] = [sem['numerator']]
    tr = pd.DataFrame()
    tr['mean'] = [ave['trace']]
    tr['standard error'] = [sem['trace']]
    covAB = cov.loc['numerator', 'trace']

    energy = ratio(num, tr, covAB, N)
    print()
    print('Energy estimate from the averaged RDM profile:')
    print(energy.to_string(float_format='%16.12E'))


def main(arguments):
    ''' Run data analysis to produce an estimate of the energy from
    the RDM numerator and RDM denominator.

    Parameters
    ----------
    arguments : list of strings
        command-line arguments.

    Returns
    -------
    None.
    '''

    filenames, options = parse_arguments(arguments)

    int_file = options.int_file
    bcomplex = options.complex
    cas = options.cas
    det = options.det

    if options.metadata_file:
        meta_data, _ = extract_data(options.metadata_file)[0]

        int_file = meta_data['system']['read_in']['int_file']
        bcomplex = meta_data['system']['read_in']['complex']
        cas = meta_data['system']['read_in']['CAS']
        det = meta_data['reference']['hilbert_space_det']

    if '/' in filenames[0]:
        path = ''
        for s in filenames[0].split('/')[:-1]:
            path += s + '/'
        int_file = path + int_file

    integrals, nel, norb, uhf = read_dump_file(int_file)

    if cas[0] > 0:
        N = cas[1]
        Mfrozen = nel - N
    else:
        N = nel
        Mfrozen = 0

    gamma = {'ijab': [], 'gamma': []}
    estimates = {'numerator': [], 'trace': []}

    for rdm in filenames:

        tmp_gamma = read_density_matrix_file(rdm)

        numerator, trace = calculate_density_matrix_energy(
                                    tmp_gamma, integrals, N, norb,
                                    Mfrozen, det, bcomplex,
                                )

        estimates['numerator'].append(numerator)
        estimates['trace'].append(trace)

        for ijab, gamma_ijab in tmp_gamma.items():
            i, j, a, b = ijab
            gamma['ijab'].append(f'{i} {j} {a} {b}')
            gamma['gamma'].append(gamma_ijab)

    stdout_average_gamma(pd.DataFrame(gamma))
    stdout_average_energy(pd.DataFrame(estimates))


if __name__ == '__main__':
    main(sys.argv[1:])
