#! /bin/bash

if [ $# -eq 0 ]; then
    cat <<END
$0: analyse a set of DMQMC RDM files.

No arguments supplied.

Usage:

$0 STEM RDM_1 RDM_2 ... RDM_N

where RDM_i is a stochastic RDM file and STEM is a file stem.
The following files are created:

STEM.rdm: collated RDM elements (ie each line contains the values for that RDM
          element given in all RDM_i files).
STEM.average_rdm: RDM matrix averaged over the individual stochastic RDMs.
STEM.entropy_estimates: estimates of entropy measures based on the RDM.
STEM.entropy_estimates: estimates of entropy measures based on the RDM.
END
    exit 1
fi

STEM=$1
shift
RDM_FILES=$@

tools=$(dirname $0)

$tools/join_rdms.sh $RDM_FILES > $STEM.rdm
$tools/average_rdm.py $STEM.rdm > $STEM.average_rdm
$tools/../dmqmc_rdm_noise/bin/rdm_entropy_noise.x $STEM.average_rdm > $STEM.entropy_estimates
$tools/rdm_entropy_stats.py $STEM.entropy_estimates > $STEM.final_estimates
