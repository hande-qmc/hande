Some geography...
=================

Files are organised in the HANDE repository as follows:

``./``
    Root directory of the program.
``bin/``
    Directory containing the compiled program, hande.x.  Created during
    compilation.
``config/``
    Directory containing the configuration input files used to generate makefiles.
``dest/``
    Directory containing the compiled object files and dependency files.  Created
    during compilation.
``documentation/``
    Directory containing documentation on the HANDE program.  The
    documentation is written in reStructured Text and can be converted
    into a wide range of output formats.
``src/``
    Directory containing the main source files.
``lib/``
    Directory containing "library" source files.  These are procedures which are
    not specific to the HANDE code but are generally useful.  Some are written
    by the authors, some are freely available (as noted in the source files).
``tools/``
    Directory containing scripts and tools for compiling, running and analysing
    output from HANDE.
``test_suite/``
    Directory containing a set of tests which HANDE should agree with.
