#!/bin/sh
#PBS -l walltime=48:00:00
#PBS -lselect=4:ncpus=12:mpiprocs=48

watchdog="/home/nb408/hubbard_fciqmc/tools/send_softexit.py 47:00:00"
application=/home/nb408/hubbard_fciqmc/bin/hubbard.cx1.optimised.x
input=heisenberg_12_neel.inp
output=$(basename $input .inp).out

workdir=$PBS_O_WORKDIR

scratchdir=$TMPDIR

module load intel-suite/11.1 mpi

cp $workdir/$input $workdir/restart.* $scratchdir

cd $scratchdir

echo "Nodes:"
cat $PBS_NODEFILE

echo "Contents of scratchdir $scratchdir:"
ls $scratchdir

echo "Running..."

$watchdog &
wps=$!
mpiexec $application $input > $output
kill -2 $wps

echo "Finished calculation."

cp $scratchdir/* $workdir

echo "All done."
