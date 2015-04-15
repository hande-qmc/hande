Adding a new test
-----------------

#.  Ensure the test suite passes with the master on your system.
#.  Now checkout the branch you're working on where you'd like to add the test.
#.  Rebuild HANDE so that the HANDE binary prints out the SHA1 hash of the current
    commit.  Make sure that there are no uncommitted changes to the source directory so
    that the benchmarks can be reproduced at a later date using the same binary.
#.  Inside test_suite find the appropriate directory in which to add your test, or
    create a new directory, appropriately named, if necessary.
#.  Inside this directory create a new directory with a sensible name describing your
    test, and change to it.
#.  Place the input files for your test in the directory.  You can have multiple input
    files in a single directory.
#.  git add your directory (this avoids having to separate out files generated during
    the tests).
#.  If you created a directory for a new category of tests then you will probably
    need to add the directory name in [ ] to the jobconfig file. If not, then the
    test should already be included through the globbing in jobconfig.
#.  If required, pick some appropriate categories to add your test to in jobconfig.
#.  Run testcode.py make-benchmarks to create new benchmarks e.g.

    .. code-block:: bash

        $ ../../testcode2/bin/testcode.py make-benchmarks
        Using executable: /home/Alex/code/HANDE/master/test_suite/../bin/hande.x.
        Test id: 09042014-2.
        Benchmark: 288ad50.

        ...

        Failed tests in:
            /home/Alex/code/HANDE/master/test_suite/H2-RHF-cc-pVTZ-Lz
        Not all tests passed.
        Create new benchmarks? [y/n] y
        Setting new benchmark in userconfig to be 6d161d0.

    Hopefully the only failed tests are your new tests (which you've checked).

    Alternatively, a better method is to make a benchmark for the new test only:

    .. code-block:: bash

        $ ../../testcode2/bin/testcode.py make-benchmarks -ic fciqmc/H2-RHF-cc-pVTZ-Lz

        ...

        Setting new benchmark in userconfig to be: 6d161d0 288ad50.

    The use of the 'i' flag tells testcode2 to insert the new benchmark at the
    start of the existing list of benchmarks, as can be seen in this example.

    If you leave the 'i' flag out then it will remove all old benchmarks, which
    we do not want.
#.  Now remember to add the benchmark files and the jobconfig and userconfig files
    to the repository.

    .. code-block:: bash

        $ git add userconfig jobconfig fciqmc/*/benchmark.out.6d161d0.inp*

    where 6d161d0 is the hash of the newly-created benchmark.

#.  Do a quick git status to make sure you haven't missed anything important out, and
    then you're ready to commit the tests:

    .. code-block:: bash

        $ git commit -m "Added new test H2-RHF-cc-pVTZ-Lz and benchmark 6d161d0."

    Remember you're committing to a branch not the master.
#.  Push this to the main repository and send round a pull request for review before its
    to be merged with master.
