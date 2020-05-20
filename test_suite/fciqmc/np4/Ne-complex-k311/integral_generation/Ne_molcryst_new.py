import numpy as np
import pyscf.pbc
from pyscf.pbc import cc as pbccc
import pyscf.tools as tools
import pyscf.pbc.tools.pbc as tools
import ase.dft.kpoints
import pyscfdump
import shutil

 
def run_kccsd(mf):
    cc = pbccc.KCCSD(mf)
    cc.verbose = 7
    cc.ccsd()
    return cc
 
def run_krccsd(mf):
    cc = pbccc.KRCCSD(mf)
    cc.verbose = 7
    cc.ccsd()
    return cc
 
def run_ip_krccsd(cc, nroots=9):
    e,c = cc.ipccsd(nroots)
    return e,c
 
def run_ea_krccsd(cc, nroots=9):
    e,c = cc.eaccsd(nroots)
    return e,c
if __name__ == '__main__':
    import sys
    from pyscfdump.helpers import get_ase_diamond_primitive, build_cell, get_ase_atom
    from ase.build import bulk
    from pyscfdump.scf import run_khf
    from pyscfdump import pbcfcidump
    import ase
 
    args = sys.argv[1:]
    if len(args) != 13:
        print '''usage: atom basis ke nkx nky nkz gamma lattice startlen endlen skiplen doqmc exxdiv
atom is the symbol of the atom to make the crystal from
basis is e.g. STO-3G
ke is the max KE cutoff (e.g. 20)
nkx nky nkz is the kpoint grid (e.g. 2 1 1 )
gamma is best 0, but if 1 then the kpoint mesh is shifted to include the gamma point explicitly (required for HANDE).
lattice is lattice-parameter a in bohr
startlen is the initial bond length for the scan (truncated to 2dp) in bohr
endlen is the final bond length for the scan (truncated to 2dp) in bohr
skiplen is the skip interval for the scan (truncated to 2dp) in bohr
doqmc:
   0 just CC
   1 CCMC & CC
   2 CCMC
   3 just save the fcidumps
exxdiv 0 = noG0 1= ws_cut
'''
        sys.exit(1)
    atom = args[0]
    bas = args[1]
    if '/' in bas:
        bas, pp = bas.split('/')
    else:
        pp = None
    ke = float(args[2])
    nmp = np.array([int(nk) for nk in args[3:6]])
    gamma= bool(int(args[6]))
    lattice=float(args[7])
    startlen=int(float(args[8])*100)
    endlen=int(float(args[9])*100)
    skiplen=int(float(args[10])*100)
    doqmc=int(args[11])
    exd=int(args[12])
    if exd==0:
      exxdiv=None
    elif exd==1:
      exxdiv='vcut_ws'
   
#    assert atom in ['C','Si']

    ds=[]
    es=[]
    cs=[]
    hes=[]
#0.7122
    import subprocess,os
    my_env = os.environ.copy()
    if "PBS_NUM_PPN" in my_env:
      my_env["OMP_NUM_THREADS"]=my_env["PBS_NUM_PPN"]

#    for i in range(135,136):

#range(60,550,5):
#    for i in [135]: #,500]:
    ks=str(nmp[0])+str(nmp[1])+str(nmp[2])
    foutp=open("ne_"+ks+"_k"+str(ke)+".ccmc.out",'w')
    for i in range(startlen,endlen,skiplen):
        lattice = i/100.0
        dist = lattice
#        ase_atom = ase.Atoms("HH",positions=[[0,0,0],[dist,0,0]],cell=[lattice,lattice,lattice],pbc=True)
#        print dist

    #    ase_atom = bulk(atom,a=lattice,crystalstructure='rocksalt') #get_ase_atom("lih") #get_ase_diamond_primitive(atom=atom) 
        ase_atom = bulk(atom,a=lattice,crystalstructure='fcc') #get_ase_atom("lih") #get_ase_diamond_primitive(atom=atom) 
        cell = build_cell(ase_atom, ke=ke, basis=bas, incore_anyway=True, pseudo=pp) #pseudo=None) 
#,crystalstructure="bcc")
     
        mf,scaled_kpts = run_khf(cell, nmp=nmp,gamma=gamma ,exxdiv=exxdiv)
        ds.append(dist)
        es.append(mf.e_tot)
        if doqmc==0 or doqmc==1:
            cc = run_krccsd(mf)
            cs.append(cc.ecc)
        if doqmc>0:
            fcid = 'fcidumpfile'
            pbcfcidump.fcidump(fcid,mf,nmp,scaled_kpts,not gamma) 
            if doqmc==3:
               shutil.copyfile(fcid,fcid+'.r'+str(dist))
               shutil.copyfile(fcid+'_X',fcid+'_X.r'+str(dist))
               d,e=zip(ds,es)[-1]
               foutp.write(str(d)+"\t"+str(e)+"\n")
               foutp.flush()
            else:
               import subprocess
               po=subprocess.Popen(["/home/ajwt3/code/HANDE/bin/hande.x","HH_ccmc.in"],stdout=subprocess.PIPE) 
               check=0
               cclist=[]
               fcc=open("nne"+ks+"_k"+str(ke)+"_r"+str(dist)+".ccmc.out",'w')
               for x in iter(po.stdout.readline, ""):
                   fcc.write(x)
                   fcc.flush()
                   if 'Correlation energy:' in x:
                       check=1
                   if check:
                       check+=1
                   if check>0 and check<6:
                       sp=x.split()
                       cclist.append(sp[-1])
               hes.append(cclist)
               fcc.close()
               d,e,c,he=zip(ds,es,cs,hes)[-1]
               foutp.write(str(d)+"\t"+str(e)+"\t"+str(c)+"\t"+str(he[0])+"\t"+str(he[1])+"\t"+str(he[2])+"\t"+str(he[3])+"\n")
               foutp.flush()
        else:
            d,e,c=zip(ds,es,cs)[-1]
            foutp.write(str(d)+"\t"+str(e)+"\t"+str(c)+"\n")
            foutp.flush()
    if doqmc==2:
        for d,e,he in zip(ds,es,hes):
             print d,e,he[0],he[1],he[2],he[3]
    elif doqmc==1:
        for d,e,c,he in zip(ds,es,cs,hes):
             print d,e,c,he[0],he[1],he[2],he[3]
    elif doqmc==3:
        for d,e in zip(ds,es):
             print d,e
    else:
        for d,e,c in zip(ds,es,cs):
             print d,e,c
