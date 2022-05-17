#!/usr/bin/env python
''' rdm_energy_estimate.py [options] RDM_1 RDM_2 ... RDM_N

Average together RDM files produced during an FCIQMC calculation and
perform additional analysis such as calculating the energy estimate from 
the RDM profile.
'''

import sys
import pkgutil
import argparse
import numpy as np

try:
    from pyhande.extract import extract_data
except ModuleNotFoundError:
    _script_dir = os.path.dirname(os.path.abspath(__file__))
    if not pkgutil.find_loader('pyhande'):
        sys.path.append(os.path.join(_script_dir, '../pyhande'))
    from pyhande.extract import extract_data


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

    parser.add_argument('-o', '--output', default='csv',
                        choices=['csv', 'txt'], help='Format for data table. '
                        '  Default: %(default)s.')
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
    parser.add_argument('-of', '--output_file', action='store', default=None,
                        type=str, dest='output_file', help='Define a file '
                        'name to store the results of averaging to.')
    parser.add_argument('-ff', '--float-format', action='store',
                        default=None, type=str, dest='float_format',
                        help='Format the values from the resulting analysis '
                        'before reporting, only used for the "txt" format. '
                        'I.E., %%6.4f, %%12.8E, ..., etc.')
    parser.add_argument('-tol', '--tolerance-number', action='store',
                        default=None, type=int, dest='csv_tol', help='Set the '
                        'tolerance of the reported data from analysis. For '
                        'example, 12 would round the final data report '
                        '(after analysis) to the 12th decimal place before '
                        'reporting. Only used for the "csv" format.')
    parser.add_argument('filenames', nargs='+', help='HANDE files to analyse.')
    parser.parse_args(args=None if arguments else ['--help'])

    options = parser.parse_args(arguments)

    return options.filenames, options


def read_dump_file(dump):

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
                line_data = line.split()
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

    return integrals, nel, norb, uhf


def read_density_matrix_file(rdm):

    gamma = {}

    with open(rdm, 'r') as stream:

        for line in stream:
            line_data = line.split()
            i, j, a, b = (int(k) for k in line_data[:-1])
            gamma_ijab = float(line_data[-1])
            gamma[i, j, a, b] = gamma_ijab

    return gamma


def orbital_index_from_pair_index(index_pair):

    index_one = int(1.5 + (2*index_pair-1.75)**0.5)
    index_two = index_pair - int(int((index_one-1)*(index_one-2))/2)

    if index_one > index_two:
        index_one, index_two = index_two, index_one

    return index_one, index_two


def get_two_particle(integrals, i, j, k, l, bcomplex):

    # http://vergil.chemistry.gatech.edu/notes/permsymm/permsymm.html

    if not bcomplex:
        symmetries = [
                [i, j, k, l],
                [j, i, l, k],
                [k, l, i, j],
                [l, k, j, i],
                [k, j, i, l],
                [l, i, j, k],
                [i, l, k, j],
                [j, k, l, i],
            ]
    else:
        symmetries = [
                [i, j, k, l],
                [j, i, l, k],
                [k, l, i, j],
                [l, k, j, i],
            ]

    for ii, jj, kk, ll in symmetries:
        if (ii, jj, kk, ll) in integrals.keys():
            return integrals[ii, jj, kk, ll]

    return 0.0


def search_twofold(integrals, i, a):

    if (i, 0, a, 0) in integrals.keys():
        return integrals[i, 0, a, 0]
    elif (a, 0, i, 0) in integrals.keys():
        return integrals[a, 0, i, 0]

    return 0.0


def get_one_particle(integrals, a, b, Nfrozen, bcomplex):

    hab = search_twofold(integrals, a, b)

    if Nfrozen == 0:
        return hab

    for i in range(1, Nfrozen + 1):
        hab += get_two_particle(integrals, i, a, i, b, bcomplex)
        hab -= get_two_particle(integrals, i, a, b, i, bcomplex)

    return hab


def calculate_density_matrix_energy(gamma, integrals, nel, norb, Mfrozen,
                                    det, bcomplex, error=False):

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

    gamma_all = {}

    for rdm in filenames:
        integrals, nel, norb, uhf = read_dump_file(int_file)

        gamma = read_density_matrix_file(rdm)

        for k, v in gamma.items():
            if k not in gamma_all.keys():
                gamma_all[k] = []

            gamma_all[k].append(v)

    if cas[0] > 0:
        N = cas[1]
        Mfrozen = nel - N
    else:
        N = nel
        Mfrozen = 0

    gamma_ave = {}
    gamma_err = {}

    for k, v in gamma_all.items():
        gamma_ave[k] = np.mean(v)
        gamma_err[k] = (np.std(v)/np.sqrt(len(v)-1))**2.0

    numerator, trace = calculate_density_matrix_energy(gamma_ave, integrals,
                                                       N, norb, Mfrozen, det,
                                                       bcomplex)

    sem_nume, sem_tr = calculate_density_matrix_energy(gamma_err, integrals,
                                                       N, norb, Mfrozen, det,
                                                       bcomplex)

    sem_nume = abs(sem_nume)

    energy = numerator/trace
    error = ((sem_nume**0.5/numerator)**2 + (sem_tr**0.5/trace)**2)**0.5
    error *= energy

    print(f' E_RDM: {energy:>18.12E} +/- {error:>18.12E}')


if __name__ == '__main__':
    main(sys.argv[1:])
