#!/usr/bin/env python

# Irreps ID maps
# Dooh	->  D2h				|		Coov	-> C2v
# A1g	0		Ag	 0		|		A1		0		 A1		 0
# A2g	1		B1g	 1		|		A2		1		 A2		 1
# A1u	5		B1u	 5		|		E1x		2		 B1		 2
# A2u	4		Au	 4		|		E1y		3		 B2		 3
# E1gx	2		B2g	 2		|		E2x		10		 A1		 0
# E1gy	3		B3g	 3		|		E2y		11		 A2		 1
# E1ux	7		B3u	 7		|		E3x		12		 B1		 2
# E1uy	6		B2u	 6		|		E3y		13		 B2		 3
# E2gx	10		Ag	 0		|		E4x		20		 A1		 0
# E2gy	11		B1g	 1		|		E4y		21		 A2		 1
# E2ux	15		B1u	 5		|		E5x		22		 B1		 2
# E2uy	14		Au	 4		|		E5y		23		 B2		 3
# E3gx	12		B2g	 2		|
# E3gy	13		B3g	 3		|
# E3ux	17		B3u	 7		|
# E3uy	16		B2u	 6		|
# E4gx	20		Ag	 0		|
# E4gy	21		B1g	 1		|
# E4ux	25		B1u	 5		|
# E4uy	24		Au	 4		|
# E5gx	22		B2g	 2		|
# E5gy	23		B3g	 3		|
# E5ux	27		B3u	 7		|
# E5uy	26		B2u	 6		|

#Psi4 (Cotton) to PySCF
#D2h
#Irrep PySCF Psi4
#Ag		 0		 1
#B1g	 1		 4
#B2g	 2		 6
#B3g	 3		 7
#Au		 4		 8
#B1u	 5		 5
#B2u	 6		 3
#B3u	 7		 2

# This is a utility that outputs Lz-transformed FCIDUMP integral files from a PySCF calculation.
# Currently only supports linear D_infh molecules. The system can be specified at the bottom of this script.
# For each single-point calculation this script will output three files:
# 1. Original untransformed PySCF output (filename-pyscf.FCIDUMP)
# 2. Untransformed output with Psi4 symmetry ordering and HF eigenvalues added (filename.FCIDUMP)
# 3. Lz-transformed output (filename-lz.FCIDUMP)

# Unlike other QC programs, PySCF is uniquely capable of handling (subsets of) the infinite point groups (up to E5) and 
# it is also capable of separating the x- and y-like component of the degenerate (E) irreps, giving
# non-standard 'irreps' like E1ux for D_infh, which would be spanned by an in-phase pi bond formed of 
# two px orbitals on either ends of a homonuclear dimer. With this information we will be able to identify which pairs of MOs
# can be transformed into eigenfunctions of the Lz operator (with well-defined Ml values). 

# Example:
# We know the standard spherical harmonic coefficients: 
#			C(1,+1) \propto -(x+iy)/sqrt(2)
#			C(1,-1) \propto +(x-iy)/sqrt(2)
# which means (-1,-1|+1,+1) = 1/2 * \iint [(x+iy)*(x-iy)] 1/r12 [-(x-iy)*-(x+iy)] dr1 dr2 
#									(remember to take complex conjugates for the 1st and 3rd functions)
#							= 1/2 * [(xx+yy|xx+yy)]
#							= 1/2 * [(xx|xx) + (yy|yy) + 2(xx|yy)]
# A general function for obtaining these coefficients for arbitrary Ml is in get_lz_idx_and_coeff below.

# This script first calls PySCF to generate the untransformed FCIDUMP. We then proceed to process the ORBSYM information, 
# For an 'irrep' of En(g/u)(x/y), we give the MO a 'destination' Ml label of + (for x)/ - (for y) n, and 0 otherwise.
# After which we overwrite the ORBSYM to only reflect the inversion symmetry, with +1 for g and -1 for u.
# We can now begin the transformation of the integrals. In general, for a 2e integral we will get 
# at most 16 contributions to an Lz-transformed integral,
# if all 4 functions were previously degenerate. We simply collect all 16 of these contributory integrals 
# and their coefficients and add them up. 1e integrals are processed identically.

# PySCF uses the sensible numbering system for D2h (see above), for now we convert this to the (less sensible) 
# Cotton's ordering used in Psi4 (see above) for the untransformed FCIDUMP.

from pyscf import gto, scf, mcscf, tools

import numpy as np
import f90nml
import subprocess as sp
import sys
import os

# Global constant for the coefficients
NORM = 1/np.sqrt(2)

def orbsym_lz_transform(_orbsym):
	"""
	Process the PySCF D_infh symmetry irrep labels and return
	_orbsym with only inversion symmetry (1 for gerade and 2 for ungerade),
	_syml with all '-20' filled in, and _symlz with the correct Ml value for the orbital.
	"""
	_ci_orbsym = np.zeros_like(_orbsym)
	_syml = np.zeros_like(_orbsym)
	_syml[:] = -20
	_symlz = np.zeros_like(_orbsym)
	for idx, item in enumerate(_orbsym):
		if item in [0,1,5,4]:
			_symlz[idx] = 0
		elif item in [2,7]:
			_symlz[idx] = 1
		elif item in [3,6]:
			_symlz[idx] = -1
		elif item in [10,15]:
			_symlz[idx] = 2
		elif item in [11,14]:
			_symlz[idx] = -2
		elif item in [12,17]:
			_symlz[idx] = 3
		elif item in [13,16]:
			_symlz[idx] = -3
		elif item in [20,25]:
			_symlz[idx] = 4
		elif item in [21,24]:
			_symlz[idx] = -4
		elif item in [22,27]:
			_symlz[idx] = 4
		elif item in [23,26]:
			_symlz[idx] = -4
		else:
			raise ValueError(f'Unsupported orbsym: {_symlz[idx]}')
		_ci_orbsym[idx] = ((item%10)>3)+1
	return _ci_orbsym,_syml,_symlz

def pyscf_to_psi4(_pyscf_orbsym):
	_dict = {0:1,1:4,2:6,3:7,4:8,5:5,6:3,7:2}
	_psi4_orbsym = np.zeros_like(_pyscf_orbsym)
	for i in range(len(_pyscf_orbsym)):
		_psi4_orbsym[i] = _dict[_pyscf_orbsym[i]%10]

	return _psi4_orbsym

def pair_lz(_symlz, _eigval):
	"""
	Find pairs of degenerate orbitals with opposive Ml, and store
	the partners' indices in a list. The indices are 0-indexed.
	"""
	_lzpairs = np.zeros_like(_symlz)
	_lzpairs[:] = -1
	for i in range(len(_symlz)):
		if _symlz[i] == 0:
			_lzpairs[i] = i
		elif _lzpairs[i] > 0:
			# Already paired
			continue
		else:
			# symlz != 0 and not paired already
			for j in range(i,len(_symlz)):
				if _symlz[i] == -_symlz[j] and abs(_eigval[i] - _eigval[j])<1e-10:
					_lzpairs[i] = j
					_lzpairs[j] = i
					break
	if any(_lzpairs==-1):
		print(_lzpairs)
		raise ValueError('Orbitals not paired correctly!')
	return _lzpairs

def tri_ind(i,j):
	"""
	The triangular index of a symmetric matrix (1-indexed)
	1 2 4
	2 3 5
	4 5 6
	"""
	if i >= j:
		a = (i*(i-1))/2 + j
	else:
		a = (j*(j-1))/2 + i
	return int(a)
	
def eri_ind(i,j,k,l):
	"""
	The two-electron ERI tensor is 8-fold symmetric (note we're using the Chemists' notation)
	(ij|kl) = (ij|lk) = (kl|ij) = (kl|ji) = (ji|kl) = (ji|lk) = (lk|ij) = (lk|ji)
	So we compress the tensor into a one-dimensional list by only storing
	i >= j, k >= l, and ij >= kl, where ij is the triangular index (defined above).
	Note that this function input is 1-indexed as per Fortran, and output is also 1-indexed as the 
	zeroth index is the default (0) value
	"""
	
	if (i==0 or j==0 or k==0 or l==0):
		return 0
	return int(tri_ind(tri_ind(i,j), tri_ind(k,l)))

def get_lz_idx_and_coeff(_lz,_i,_pair,_conj):
	"""
	Find out which orbitals contribute to an Lz-transformed orbital,
	and their coefficients.
	In:
		_lz: The signed Ml value.
		_i: The (0-indexed) index of the current real orbital
		_pair: The (0-indexed) index of the Lz-paired orbital of the current real orbital/
		_conj: whether we should conjugate this Lz-transformed orbital: 1 for no, -1 for yes

	Return:
		_a1: The (0-indexed) index of the x-like (positive Ml) contributing orbital.
		_a2: The (0-indexed) index of the y-like (negative Ml) contributing orbital, zero if _lz=0.
		_a1c: The coefficient of _a1.
		_a2c: The coefficient of _a2.
	"""
	# (ij|kl) = \iint i*(1)j(1)(1/r12)k*(2)l(2)
	if _lz == 0:
		_a1 = _i
		_a2 = 0
		_a1c=1+0j
		_a2c=0+0j
	else:
		if _lz < 0:
			_a1 = _pair
			_a2 = _i
			# odd and even are the same coefficients
			_a1c = complex(NORM,0)
			_a2c = complex(0,-NORM*_conj)
		else:
			_a1 = _i
			_a2 = _pair
			if (abs(_lz)%2 == 1):
				# Odd Ml
				_a1c = complex(-NORM,0)
				_a2c = complex(0,-NORM*_conj)
			else:
				# Evem Ml
				_a1c = complex(NORM,0)
				_a2c = complex(0,NORM*_conj)	 
	return _a1, _a2, _a1c, _a2c

def run_pyscf(_mol, _occ=None):
	# Run RHF
	_mf = scf.RHF(_mol)
	if (_occ is not None):
		# Clamp occupancy to improve convergence (tracking is another option)
		_mf.irrep_nelec = _occ
	_mf.kernel()

	return _mf

def make_fcidumps(_mf, _filename):
	# Dump out FCIDUMP, format string is such that it produces the same output as Psi4
	print(f'Generating untransformed PySCF FCIDUMP in {_filename}-pyscf.FCIDUMP...')
	tools.fcidump.from_scf(_mf, f'{_filename}-pyscf.FCIDUMP', tol=1e-12, float_format='%28.20E')
	print('Done!')

	# PySCF doesn't print out HF eigenvalues, so we do it ourselves
	eigval = _mf.mo_energy
	e_nuc = _mf.energy_nuc()

	# f90nml by default reads the whole FCIDUMP in, which can be slow, 
	# here we just 'grep' the namelist bit and pass to f90nml
	nmlstr = ''
	with open(f'{_filename}-pyscf.FCIDUMP','r') as f:	
		while True:
			line = f.readline()
			nmlstr += line
			if ('&END' in line):
				break

	# Read the Fortran namelist, we can also get these variables from the mf object
	#nml = f90nml.read(f'{_filename}-pyscf.FCIDUMP')
	nml = f90nml.reads(nmlstr)
	norb = int(nml['fci']['norb'])
	nelec = int(nml['fci']['nelec'])
	pyscf_orbsym = np.array(nml['fci']['orbsym'],dtype='int')

	# We print out the Psi4 ORBSYM
	psi4_orbsym = pyscf_to_psi4(pyscf_orbsym)

	# Process orbsym to get Lz information, and then only keep inversion
	# symmetry in orbsym (1 for gerade and 2 for ungerade)
	orbsym, syml, symlz = orbsym_lz_transform(pyscf_orbsym)

	# Find out which orbitals are degenerate
	lzpairs = pair_lz(symlz, eigval)

	# Allocate arrays
	pair = int((norb*(norb+1))/2)
	size = int((pair*(pair+1))/2)
	print(f'Max size of 2e integrals: {size}')
	# zero-indexed
	oeint = np.zeros((norb,norb))
	# zeroth element is zero, so the rest of the array is Fortran style/1-indexed
	teint = np.zeros(size+1)

	# Read the FCIDUMP in, skipping the namelist
	read = False
	print(f'''Reading PySCF FCIDUMP and writing the Psi4-compatible FCIDUMP 
		(with Cotton D2h irrep ordering and HF eigenvalues) in {_filename}.FCIDUMP''')
	with open(f'{_filename}-pyscf.FCIDUMP', 'r') as infile, open(f'{_filename}.FCIDUMP', 'w') as outfile, open(f'{_filename}-temp.FCIDUMP', 'w') as tempfile:
		# Write tempfile namelist

		# Compatible format with Psi4
		tempfile.write(f'&FCI\nNORB={int(norb)},\nNELEC={int(nelec)},\nMS2=0,\nORBSYM=')
		for i in orbsym:
			tempfile.write(f'{int(i)},')
		tempfile.write('\n')
		tempfile.write('ISYM=1,\nUHF=.FALSE.,\nSYML=')
		for i in syml:
			tempfile.write(f'{int(i)},')
		tempfile.write('\nSYMLZ=')
		for i in symlz:
			tempfile.write(f'{int(i)},')
		tempfile.write('\nPAIR=')
		for i in lzpairs:
			# lzpairs are 0-indexed, make it Fortran-compatible
			tempfile.write(f'{int(i+1)},')
		tempfile.write('\n&END\n')

		contents = infile.readlines()
		for line in contents:
			if read:
				data = line.strip().split()
				i = int(data[1])
				j = int(data[2])
				k = int(data[3])
				l = int(data[4])

				# Add HF eigenvalues as PySCF doesn't support writing them out
				if i==0 and j==0 and k==0 and l==0:
					for i in range(norb):
						outfile.write(f'{eigval[i]:28.20E}{(i+1):4d}{0:4d}{0:4d}{0:4d}\n')
						tempfile.write(f'{eigval[i]:28.20E}{(i+1):4d}{0:4d}{0:4d}{0:4d}\n')
					# Don't forget E_core
					outfile.write(line)
					tempfile.write(line)
					break
				else:
					tempfile.write(line)

			if '&END' in line:
				read = True

			# Write exactly as read in, but translate ORBSYM to Psi4 order
			if 'ORBSYM' in line:
				outfile.write('ORBSYM=')
				for i in psi4_orbsym:
					outfile.write(f'{int(i)},')
				outfile.write('\n')
			else:
				outfile.write(line)

	print('Handing over to Fortran for the Lz transform...')

	sp.run(['./lz_transform.x', f'{_filename}-temp.FCIDUMP'], stderr=sys.stderr, stdout=sys.stdout)
	sp.call(f'rm {_filename}-temp.FCIDUMP', shell=True)
	sp.call(f'mv lz-{_filename}-temp.FCIDUMP {_filename}-lz.FCIDUMP', shell=True)

	print('Done!')

if __name__ == '__main__':

	if not os.path.exists('./lz_transform.x'):
		raise OSError('lz_transform.x does not exist, please compile lz_transform.f90!')

	#for bl in np.linspace(2,7,11):
	for bl in [2.5]:

		# PySCF FCIDUMP name
		filename = f'be-ccpvdz-{bl:.1f}'

		print('='*30+'\n'+f'Doing {filename}\n'+'='*30)

		# Define molecule
		mol = gto.M(atom=f"""
		Be
		Be 1 {bl:.1f}""",
		basis='cc-pvdz', symmetry='Dooh')

		# Clamp occupancy for better convergence / converging to the same RHF solution across a binding curve
		# Restarting from previous geometry (ie tracking) works too.
		occ = {'A1g':4,'A1u':4}

		mf = run_pyscf(mol, occ)

		make_fcidumps(mf, filename)