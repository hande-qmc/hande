#!/bin/bash

usage() {
    cat <<END
convert_input.sh file1 file2 ... fileN

Convert set of input files to use the "new" input options for specifying the system type:

k_space -> hubbard_k
momentum_space -> hubbard_momentum (identical to hubbard_k)
real_space -> hubbard_real
END
}

if [[ $# -eq 0 ]]; then
    usage
    exit 1
fi

for input in $@; do

    # New keywords for system type?  If so, do nothing.
    if grep -qi "hubbard_k\|hubbard_momentum\|hubbard_real\|heisenberg" $input; then

        echo "New-style input file.  Not changing."
        
    # Any system type set?  If not, set hubbard_k to be the system type (explicit
    # is better than implicit!)
    elif ! grep -qi "real_space\|k_space\|momentum_space" $input; then

        echo "Setting system_type to explicitly be hubbard_k in $input."

        input_contents=$(cat $input)

        echo -e "hubbard_k\n\n$input_contents" > $input

    # Update keywords
    else

        echo "Changing $input to use new keywords."

        perl -pi -e 's/real_space/hubbard_real/i; s/k_space/hubbard_k/i; s/momentum_space/hubbard_momemtum/i' $input

    fi

done
