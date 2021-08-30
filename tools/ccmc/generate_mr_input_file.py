#!/usr/bin/env python

import numpy as np
import string
import itertools
import argparse

def exlvl(f1,f2):
    """
    Finds the excitation level between two bitstrings of the same number of electrons.
    Essentially half the Hamming distance.
    E.g. exlvl(0b11110000,0b00001111) = 4
    """
    return int(bin(f1^f2).count("1")/2)

def k_set_bits_in_n_bits(k,n):
    """
    From https://stackoverflow.com/a/40449676
    Generates all bitstrings with n electrons in k spinorbitals.
    k_set_bits_in_n_bits(2,4) = [0b0011,0b0110,0b1100,0b1010,0b0101,0b1001]
    """
    powers = [1 << e for e in range(n)]
    return [sum(bits) for bits in itertools.combinations(powers, k)]

def interleave(alpha, beta):
    """
    interleave (1011,1101) = 11011011
    https://stackoverflow.com/a/39490836
    """
    
    # Works for a max CAS of (8,8)
    #masks = [0x5555,0x3333,0xF0F]
    #powers = [1,2,4]
    #for i in range(2,-1,-1):
    
    # Works for a max CAS of (16,16)
    #masks = [0x55555555, 0x33333333, 0x0F0F0F0F, 0x00FF00FF]
    #powers = [1,2,4,8]
    #for i in range(3,-1,-1):
    
    # Works for a max CAS of (32,32)
    masks = [0x5555555555555555, 0x3333333333333333, 0xF0F0F0F0F0F0F0F0F, 0xFF00FF00FF00FF, 0xFFFF0000FFFF]
    powers = [1,2,4,8,16]
    for i in range(4,-1,-1):
        alpha = (alpha | (alpha << powers[i])) & masks[i]
        beta = (beta | (beta << powers[i])) & masks[i]
    return beta | (alpha << 1)

def ref_string(bstring,nfrozen):
    string = '{'
    for i in range(1,nfrozen+1):
        string += f'{i},'
    bit_length = len(bin(bstring)) - 2
    for shift in range(bit_length):
        if ((bstring >> shift) & 1):
            string += f'{nfrozen+shift+1},'
    string = string[:-1]
    string += '},'
    return string

def excit_generator(exlvl,e,o):
    """
    Generates all determinants within exlvl of 
    the top and bottom determinants of an active space of (e,o).
    
    For now limited to even numbers of e and o
    
    We work with an alpha string and a beta string and interleave them together in the end.
    
    E.g. excit_generator(2,8,8) will create all doubles coming from
    (0000000011111111) and (1111111100000000)
    """
    nvirt = int((o*2-e)/2) # spatial! this is alpha/beta string
    nocc = int(e/2) # ditto
    bottomstring = int(2**(e/2)-1)
    topref = int(2**e-1)<<int(o*2-e)
    def pure_excit(lvl: int):
        """
        The following inner functions generate S,D,T,.. for an alpha/beta string, so we don't need to consider combinatorics
        e.g. an overall doubles can be both alpha, both beta, or one alpha one beta.
        """
        virt = k_set_bits_in_n_bits(lvl,nvirt)
        # we're assuming the CAS is not necessarily symmetrical, otherwise can just spin flip the above list
        occ = k_set_bits_in_n_bits(nocc-lvl,nocc)
        excit_list = []
        for v,o in itertools.product(virt,occ):
            excit_list.append((v<<nocc)+o)
        return excit_list
    
    pure_singles = pure_excit(1)
    singles_return = []
    for single in pure_singles:
        singles_return.append(interleave(bottomstring,single)) # beta -> beta
        singles_return.append(interleave(single,bottomstring)) # alpha -> alpha
    top = []
    for i in singles_return:
        top.append(bit_reflect(int(o*2),i))
    singles_return += top
    if exlvl == 1:
        singles_return.append(topref)
        return singles_return
    
    pure_doubles = pure_excit(2)
    
    doubles_return = []
    for a,b in itertools.product(pure_singles,repeat=2):
        doubles_return.append(interleave(a,b))
    for double in pure_doubles:
        doubles_return.append(interleave(bottomstring,double)) # 2beta -> 2beta
        doubles_return.append(interleave(double,bottomstring)) # 2alpha -> 2alpha    
    top = []
    for i in doubles_return:
        top.append(bit_reflect(int(o*2),i))
    
    doubles_return += singles_return + top
    if exlvl == 2:
        doubles_return.append(topref)
        return doubles_return
    
    pure_triples = pure_excit(3)
    
    triples_return = []
    for a,b in itertools.product(pure_singles,pure_doubles):
        triples_return.append(interleave(a,b)) # abb -> abb
        triples_return.append(interleave(b,a)) # aab -> aab
    for triple in pure_triples:
        triples_return.append(interleave(bottomstring,triple)) # aaa -> aaa
        triples_return.append(interleave(triple,bottomstring)) # bbb -> bbb
    top = []
    for i in triples_return:
        top.append(bit_reflect(int(o*2),i))
    
    triples_return += doubles_return + top
    if exlvl == 3:
        triples_return.append(topref)
        return triples_return
    
    pure_quads = pure_excit(4)
    
    quads_return = []
    for a,b in itertools.product(pure_doubles,pure_doubles):
        triples_return.append(interleave(a,b)) # aabb -> aabb
    for a,b in itertools.product(pure_singles,pure_triples):
        triples_return.append(interleave(a,b)) # abbb -> abbb
        triples_return.append(interleave(b,a)) # aaab -> aaab
    for quad in pure_quads:
        quads_return.append(interleave(bottomstring,quad)) # aaaa -> aaaa
        quads_return.append(interleave(quad,bottomstring)) # bbbb -> bbbb
        
    top = []
    for i in quads_return:
        top.append(bit_reflect(int(o*2),i))
    
    quads_return += triples_return + top
    if exlvl == 4:
        quads_return.append(topref)
        return quads_return

def bit_reflect(blength,bstring):
    """
    bit_reflect(8,00100111)=11100100
    """
    r_bstring = 0
    for i in range(blength):
        r_bstring+= ((bstring >> i) & 1) << (blength - 1 - i)
    return r_bstring

def uncompressed(e,o,nfrozen,exlvl):
    n = int(o)
    k = int(e/2)
    bit_8 = k_set_bits_in_n_bits(k,n)
    ref_list = list(map(interleave,itertools.product(bit_8,repeat=2)))[1:]
    return ref_list

def compressed(e,o,nfrozen,exlvl):
    return excit_generator(exlvl,e,o)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generates mrccmc HANDE input files.')
    parser.add_argument("-i","--intdump", dest="intdumpname")
    parser.add_argument("-e","--nel", dest='e', help='number of electrons')
    parser.add_argument("-r","--norb", dest='o', help='number of spatial orbitals')
    parser.add_argument("-o", dest='inputname')
    parser.add_argument("--compress",default=True,type=bool,dest='compress')
    parser.add_argument("-f","--nfrozen",dest='nfrozen')
    parser.add_argument("--secref_file",dest='secref_file',help='Name of separate secondary reference file')
    parser.add_argument("--mr_excit_lvl",dest='mr_excit_lvl',help='Allowable excitation level from each reference.')
    parser.add_argument("-l","--cclvl",dest='cclvl',help='Reference CC level, 3 for CCSDT(N2), 4 for CCSDTQ(C2), 6 for CCSDTQ56(Cr2). Should be half the number of CAS orbitals.')
    args = parser.parse_args()
    intdumpname = args.intdumpname
    inputname = args.inputname
    e = int(args.e)
    o = int(args.o)
    genlvl = int(args.cclvl)-2
    compress = args.compress
    mr_n_frozen = int(args.nfrozen)
    secref_file = args.secref_file
    mr_excit_lvl = int(args.mr_excit_lvl)
    if compress:
        ref_gen = compressed
    else:
        ref_gen = uncompressed

    template = string.Template(
"""sys = read_in {
    int_file = "${intdumpname}",
    sym = tot_sym,
    Lz = true,
}
 
ccmc {
    sys = sys,
    qmc = {
        tau = 0.005,
        real_amplitudes = true,
        init_pop = 2000,
        mc_cycles = 2,
        nreports = 1e5,
        target_population = 1e5,
        state_size = -4000,
        spawned_state_size = -4000,
        vary_shift = false,
        vary_shift_from = 'proje',
        shift_damping = 0.01
        },
    restart = {
        write = 9,
        },
    reference = {
        ex_level = 2,
        },
    ccmc = {
        full_non_composite = true,
        even_selection = true,
        multiref = true,
        mr_read_in = true,
        mr_secref_file = "${mr_secref_file}",
        mr_excit_lvl = ${mr_excit_lvl},
        mr_n_frozen = ${mr_n_frozen},
        n_secondary_ref = ${nref},
     },
}
    """)

    ref_list = ref_gen(e,o,mr_n_frozen,genlvl)

    with open(secref_file,'w') as f:
        for item in ref_list:
            f.write(f'{item}\n')

    d = {'nref': len(ref_list),'mr_secref_file': secref_file ,'mr_excit_lvl': mr_excit_lvl, 'mr_n_frozen':mr_n_frozen, 'intdumpname':intdumpname}
    a = template.substitute(d)
    with open(inputname,'w') as f:
        f.write(a)
