#!/usr/bin/python

# WARNING: this is a hack.  Only common (simple) cases are considered.

import sys
import collections

def read_old(filename):
    '''Read the old input file in and convert it to a dictionary.'''

    inp = {}
    lattice = False
    idim = 0
    ndim = 3
    with open(filename) as f:
        for line in f:
            words = line.lower().split()
            if not words:
                pass
            elif words[0] in ('read', 'dipole_integrals'):
                # don't lower-case filenames!
                inp[words[0]] = line.split()[1:]
            elif words[0] == 'lattice':
                inp[words[0]] = []
                lattice = True
            elif lattice and idim != ndim :
                ndim = len(words)
                idim += 1
                inp['lattice'].append(words)
            elif len(words) == 1:
                inp[words[0]] = None
            else:
                inp[words[0]] = words[1:]
    return inp

def read_to_dict(inp, keys, remap, store):
    '''Select keys (and possibly remap them) from the input dictionary and copy them to the store.'''

    for key_orig in keys:
        # rename keys
        if key_orig in remap.keys():
            key = remap[key_orig]
        else:
            key = key_orig
        if key_orig in inp.keys():
            # hard-coded special cases
            if key_orig == '2d':
                store['dim'] = [2]
            elif key_orig == '3d':
                store['dim'] = [3]
            # general
            elif inp[key_orig]:
                store[key] = inp[key_orig]
            else:
                store[key] = ['true']
    return store

def get_sys(inp):
    '''Get the system settings from the input dictionary.'''

    system = collections.OrderedDict()

    for sys in 'hubbard_real hubbard_k heisenberg read ueg chung-landau'.split():
        if sys in inp.keys():
            system['type'] = sys
            break
    if system['type'] == 'chung-landau':
        system['type'] = 'chung_landau'
    if system['type'] == 'read':
        system['type'] = 'read_in'
        if inp['read']:
            system['int_file'] = inp['read']

    system_keys = ('electrons nel lattice ms sym U t ktwist finite_cluster '
                   'triangular_lattice Lz cas 2d 3d ecutoff rs chem_pot '
                   'dipole_integrals'.split())

    remap = {
                'finite_cluster':'finite',
                'ecutoff':'cutoff',
                'dipole_integrals':'dipole_int_file',
                'twist':'ktwist',
            }

    read_to_dict(inp, system_keys, remap, system)

    return system

def get_calc(inp):
    '''Get the calculation settings from the input dictionary.'''

    calc = collections.OrderedDict()

    remap = {
                'estimate_canonical_kinetic_energy':'kinetic_energy',
                'estimate_hilbert_space':'hilbert_space',
            }

    for calc_type in 'estimate_hilbert_space estimate_canonical_kinetic_energy'.split():
        if calc_type in inp.keys():
            if calc_type in  remap.keys():
                calc['type'] = remap[calc_type]
            else:
                calc['type'] = calc_type
            break

    if calc['type'] == 'hilbert_space':
        calc['ncycles'] = inp['estimate_hilbert_space']

    remap = {
                'truncation_level':'ex_level',
                'seed':'rng_seed',
                'num_kinetic_cycles':'ncycles',
            }

    calc_keys = 'truncation_level seed reference num_kinetic_cycles fermi_temperature init_beta init_pop'.split()

    read_to_dict(inp, calc_keys, remap, calc)

    if calc['type'] == 'kinetic_energy':
        for (old, new) in (('init_beta', 'beta'), ('init_pop', 'nattempts')):
            if old in calc.keys():
                calc[new] = calc[old]
                calc.pop(old)

    return calc

def dict_to_table(d):
    '''Convert a python dictionary to a lua table in string format, minus the {} delimiters.'''

    s = []
    for (key, val) in d.items():
        if len(val) == 1:
            try:
                x = float(val[0])
                s.append('    %s = %s,' % (key, val[0]))
            except ValueError:
                if val[0] in ['true', 'false']:
                    s.append('    %s = %s,' % (key, val[0]))
                else:
                    s.append('    %s = "%s",' % (key, val[0]))
        elif key == 'lattice':
            lattice = '{ {%s} }' % ('}, {'.join(', '.join(x) for x in val))
            s.append('    %s = %s,' % (key, lattice))
        else:
            # have a (numerical) vector.
            s.append('    %s = {%s},' % (key, ', '.join(val)))
    return '\n'.join(s)

def print_new(sys, calc):
    '''Print out the input file in lua format.'''

    print('sys = %s {' % sys['type'])
    sys.pop('type')
    print(dict_to_table(sys))
    print('}\n')

    print('%s {' % calc['type'])
    calc.pop('type')
    print('    sys = sys,')
    print(dict_to_table(calc))
    print('}')

if __name__ == '__main__':

    inp = read_old(sys.argv[1])
    sys = get_sys(inp)
    calc = get_calc(inp)
    print_new(sys, calc)
