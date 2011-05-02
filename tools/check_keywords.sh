#!/bin/bash

# A simple scipt to check that all the input options are tested and documented.

# directory location of this script relative to working directory.
script_location=$(dirname $0)

# hence location of src directory.
src=$script_location/../src

# hence location of the documentation directory.
docs=$script_location/../documentation

# hence location of the test_suite directory.
test_suite=$script_location/../test_suite

parse_source_files=$src/parse_input.F90

# keywords are found in the input parser in the format:
#  case('KEYWORD')
# and
#  case('KEYWORD1','KEYWORD2')

# This perl script returns the list of keywords.
all_keywords=$(perl -nle \
    "if (/^ *case\(\".*\"\)|^ *case\('.*'\)/i) { # select lines which contain the keywords.
        s/ //g;                                  # remove whitespace.
        s/!.*$//;                                # remove comments.
        s/case\(|\)|'|\"//gi;                    # remove the case statement and brackets.
        s/,/\n/g;                                # split keywords occuring in the same case statement into separate lines.
        print;                                   # print "naked" keyword(s)
    }" \
$parse_source_files)

# This perl script returns the list of keywords with unique meanings as regexes.
unique_keywords=$(perl -nle \
    "if (/^ *case\(\".*\"\)|^ *case\('.*'\)/i) { # select lines which contain the keywords.
        s/ //g;                                  # remove whitespace.
        s/!.*$//;                                # remove comments.
        s/case\(|\)|'|\"//gi;                    # remove the case statement and brackets.
        s/,/\\\|/g;                              # keywords occuring in the same case statement are joined to form a single regex.
        print;                                   # print "naked" keyword(s)
    }" \
$parse_source_files)

# test to see if keyword is mentioned in the documentation.
if [ -e $docs ]; then

    doc_files=$(find $docs -name \*.rst)

    if [ -z "$doc_files" ]; then
        echo "Documentation files cannot be found."
    else
        undocumented_keywords=''
        for keyword in $all_keywords; do
            undocumented_keywords="$undocumented_keywords $(grep -qwi "\*\*$keyword\*\*" $doc_files 2> /dev/null || echo $keyword)"
        done

        # Remove trailing blanks
        undocumented_keywords=$(echo $undocumented_keywords | sed -e 's/ *$//')

        # Output
        if [ ! -z "$undocumented_keywords" ]; then
            echo -e "Undocumented keywords are:\n"
            for keyword in $undocumented_keywords; do
                echo $keyword
            done
            echo
        else
            echo "All keywords appear in the documentation."
        fi
    fi

else

    echo "Documentation directory does not exist."

fi

# test to see if keyword (or a synonym) is used in the test suite.
if [ -e $test_suite ]; then

    inp_files=$(find $test_suite -name \*.inp)

    if [ -z "$inp_files" ]; then
        echo "Test suite input files cannot be found."
    else
        untested_keywords=''
        for keyword in $unique_keywords; do
            untested_keywords="$untested_keywords $(grep -qwi $keyword $inp_files 2> /dev/null || echo $keyword)"
        done

        # Remove trailing blanks
        untested_keywords=$(echo $untested_keywords | sed -e 's/ *$//')

        # Output
        if [ ! -z "$untested_keywords" ]; then
            echo -e "Untested keywords are:\n"
            for keyword in $untested_keywords; do
                echo $keyword | sed -e 's/\\|/ or /g'
            done
            echo
        else
            echo "All keywords (or their synonyms) appear in the input files in the test suite."
        fi
    fi

else

    echo "Test suite directory does not exist."

fi
