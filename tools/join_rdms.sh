#!/bin/bash

if [ $# -eq 0 ]; then
    cat <<END
$0: join together an arbitrary number of DMQMC RDM files.
All estimates of a given matrix element are placed on the same line in the
output.

No arguments supplied.

Usage:

$0 RDM_1 RDM_2 ... RDM_N > RDM

where RDM_i is a stochastic RDM file.
END
    exit 1
fi

rdmtmp=$(mktemp)
worktmp=$(mktemp)

for file in $@; do
    if [ -s $rdmtmp ]; then
        join $rdmtmp $file > $worktmp
        mv $worktmp $rdmtmp
    else
        cp $file $rdmtmp
    fi
done

cat $rdmtmp

rm $rdmtmp $worktmp
