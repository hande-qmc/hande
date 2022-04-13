#!/usr/bin/env python

# This script generates FCIDUMP files with Lz symmetry enabled.
# Your system must belong to the D_infh point group for Lz to be a good quantum number.

from pyscf import gto, scf, tools, fci, mp, cc, mcscf

#def make_fcidump(mol, r):
for r in [2.0]:
	mol = gto.Mole()
	mol.atom = [['C', (0, 0, 0)], ['C', (0, 0, r)]]
	mol.basis = 'cc-pvdz'
	# Default unit is Angstrom
	# mol.unit = 'Bohr'
	mol.symmetry = True # if you want to take molecular symmetry into account
	# Force descent into D2h group
	# mol.symmetry_subgroup= 'D2h'
	mol.build()
	
	if mol.groupname != 'Dooh':
		raise ValueError(f'The system has point group symmetry {mol.groupname}, not D_{{\\inf h\\}}')
	
	mf = scf.RHF(mol) # could use UHF/ROHF
	mf.max_cycle = 200
	mf.kernel() # runs HF
	tools.fcidump.from_mo(mol, 'FCIDUMP.C2.cc-pvdz-'+str(r)+'A', mf.mo_coeff) # dumps out integrals to the filename specified

