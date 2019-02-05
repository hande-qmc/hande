from ase.lattice import bulk
import pyscfdump.scf as scf
from pyscfdump.helpers import build_cell
from pyscfdump.pbcfcidump import fcidump

A2B = 1.889725989	#angstrom to bohr conversion

a=3.567		#lattice parameter
ke=1000		#kinetic energy cutoff in units of rydberg
basis='sto-3g'	#basis set choice
nmp = [2,1,1]	#k point sampling

#prepare the cell object
ase_atom = bulk('C', 'diamond', a=a*A2B)
cell = build_cell(ase_atom, ke=ke, basis=basis, pseudo='gth-pade')

#run the HF calculation
kmf,scaled_kpts = scf.run_khf(cell, nmp=nmp, exxdiv='ewald', gamma=True)

#dump the integrals
fcidump('fcidumpfile',kmf,nmp,scaled_kpts,False)

