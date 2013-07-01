#!/bin/bash

rm -f rdm.temp

for file in rdm.* ; do
    if [ -a rdm.temp ]; then
        join rdm.temp $file > rdm
        cp rdm rdm.temp
    else
        cp $file rdm.temp
    fi
done

rm rdm.temp
