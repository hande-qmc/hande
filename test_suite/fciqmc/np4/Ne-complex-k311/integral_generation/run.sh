#!/bin/bash
#SBATCH -C old
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 12
source ~/code/HANDE/modules
source ~/code/pyscf/loadmodules
echo "At $PWD"
export OMP_NUM_THREADS=12
export MKL_NUM_THREADS=12
if [ "$SLURM_JOBID" == "" ] ; then
   SLURM_JOBID=$( date +%Y%m%d-%H%M%S )
fi
echo "$(date) S $SLURM_JOBID $PWD" >>/home/ajwt3/calc/pyscf/JOBS
if [ "$1" == "" ] ; then

   python ../Ne_molcryst_new.py Ne gth-dzvp/gth-pade 500 3 1 1 1 10 10.5 20 .5 1 1 &> ne_311_k500.scf.out
   echo "$(date) E $SLURM_JOBID $PWD">>/home/ajwt3/calc/pyscf/JOBS
fi
