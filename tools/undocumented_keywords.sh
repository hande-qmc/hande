#!/bin/bash

# directory location of this script relative to working directory.
script_location=$(dirname $0)

# hence location of src directory.
src=$script_location/../src

parse_source_files=$src/parse_input.F90

# keywords are found in the input parser in the format:
#  case('KEYWORD')
# and
#  case('KEYWORD1','KEYWORD2')
# the perl script returns the list of keywords.
keywords=$(perl -nle \
    "if (/^ *case\('.*'\)/i) { # select lines which contain the keywords.
        s/ //g;                # remove whitespace.
        s/case\(|\)|'//g;      # remove the case statement and brackets.
        s/,/\n/g;              # split keywords occuring in the same case statement into separate lines.
        print;                 # print "naked" keyword(s)
    }" \
$parse_source_files)

# Once there is documentation, test that the keywords are present.

echo -e "Undocumented keywords are:\n"
for keyword in $keywords; do
    echo $keyword
done
