#!/bin/bash

usage () {
    cat <<END
usage: 1d_hubbard_gs.sh N U

Find the exact ground state (in terms of t) of the 1D Hubbard model at
U consisting of N sites in the simulation cell in the thermodynamic limit.
Requires Mathematica.
END
}

tools_dir=$(dirname $0)

if [[ ! $# -eq 2 ]]; then
    usage
    exit 1
fi

HUBN=$1
HUBU=$2

math_script=$(mktemp)

# so hacky
cat >$math_script <<END 
Get["$tools_dir/wu_lieb.m"]
Print[E0[$HUBN,$HUBU]]
Exit[]
END

echo -n "E[N=$HUBN,U=$HUBU] = "
math -noprompt -run "<<$math_script"

rm $math_script
