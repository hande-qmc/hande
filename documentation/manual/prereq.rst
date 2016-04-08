Prerequisites
=============

HANDE builds upon several well-written, efficient libraries to aid portability,
efficiency and sustainability.

Dependencies
------------

Fortran and C compilers
    HANDE is written in (mostly) Fortran 2003 with some C code.  We have tested HANDE
    using GCC, Intel, Cray and IBM compilers and are interested in hearing of use with
    other compilers.

    .. note::

        HANDE is relatively aggressive in adopting new language features and hence
        requires a fairly modern Fortran compiler.  In particular, gfortran 4.5 or earlier
        is unlikely to successfully compile HANDE.

LAPACK and BLAS
    Available from http://www.netlib.org/lapack/ and http://www.netlib.org/blas/ and
    vendor implementations.  Typically installed on HPC systems and available from package
    manager.  This is only required for the FCI functionality in HANDE; the performance of
    the QMC algorithms do not depend upon the quality of the LAPACK and BLAS libraries
    used.
lua 5.2
    Lua (available from http://www.lua.org) is required.  HANDE links to the lua library,
    which is used for parsing the input file.  No performance critical code is written in
    lua.

    .. note::

        The version of the AOTUS library included with HANDE is only compatible with lua
        5.2.  Later versions of AOTUS, which HANDE should also work with, support lua 5.3
        (but not 5.2 due to API changes).

MPI (parallel compilation only)
    MPI 2 is required.  We have used a variety of implementations (including OpenMPI and
    various vendor implementations).
scalapack (parallel compilation only)
    Available from http://www.netlib.org/scalapack/ and vendor implementations.  Often
    already installed on HPC systems, included in Intel Maths Kernel Library and can be
    installed from most package managers.
python 2.7+ or python 3.2+
    Almost all tools packaged with HANDE are written in python.

    .. note::

        python 2.6 or earlier python 3 versions **may** be sufficient but will probably
        require additional work.  In particular, the argparse module (included from 2.7
        and 3.2 onwards) is required and installing (especially recent versions of )
        pandas  may be problematic.  Using a recent version of python is highly
        recommended.
pandas 0.14.1+
    The HANDE data analysis tools build heavily upon the python scientific
    stack.  In particular, pandas (available from http://pandas.pydata.org) is required
    for the ``pyhande`` module and analysis scripts, almost all of which build upon
    ``pyhande``.  pandas is not required for running HANDE but is highly recommended for
    data analysis (though strictly speaking is only required if ``pyhande`` is used,
    either directly or via analysis scripts).

Bundled dependencies
--------------------

AOTUS
    AOTUS provides a nice Fortran wrapper to Lua's C-API.  For convenience (given that
    module files are Fortran-specific), AOTUS is included in the HANDE source
    distribution.

Optional dependencies
---------------------

The following are optional depedencies which add useful (in some cases almost critical)
functionality.  However, they are less likely to be compiled on HPC systems so for ease of
testing the functionality which depends upon them can be disabled at compile-time.

HDF5
    HDF5 is a library for storing scientific data and is used in HANDE for checkpointing
    (i.e. writing and reading restart files) in QMC calculations.

    Highly recommended.  Disabling HDF5 removes the ability to perform any checkpointing.

    .. note::

        HANDE requires the Fortran 2003 interface to HDF5, which is not compiled by
        default (see below), as this offers substantial advantages when working with
        dynamically sized arrays containing variables of arbitrary kinds/precision.

libuuid
    Provenance of a calculation, and the output file(s) produced by it, is an important
    topic, currently the subject of much debate in computational science.  HANDE generates
    a universally unique identifier (UUID), which is included in all files it produces.

    Highly recommended but can be disabled without impacting on performance (but perhaps
    not on the user's sanity).

    .. note::

        Some Linux distributions install libuuid but require an additional package (e.g.
        uuid-dev) to be installed in order for libuuid to exist on default search paths.
        Some luck may be found by looking under /lib or /lib64 instead of /usr/lib and
        /usr/lib64.
TRLan
    Required for (and only for) performing FCI calculations using the Lanczos algorithm.
    Available from http://crd-legacy.lbl.gov/~kewu/trlan.html.

Compilation and installation notes
----------------------------------

Some notes on compiling the less common dependencies.

lua
^^^

Lua is straightforward to compile.  For example:

.. code-block:: bash

    $ wget -O - http://www.lua.org/ftp/lua-5.2.4.tar.gz | tar xvzf -
    $ cd lua-5.2.4
    $ make linux
    $ make install INSTALL_TOP=$HOME/local

will install the lua program and library to subdirectories in $HOME/local.  It is usually
fine to compile lua using the GCC compiler and link HANDE against it using another
compiler family (e.g. Intel).

HDF5
^^^^

HDF5 uses the GNU autotools build system, so is also straightforward to compile.  For
example:

.. code-block:: bash

    $ wget -O - http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.14.tar.gz | tar xvzf -
    $ cd hdf5-1.8.14
    $ ./configure --prefix=$HOME/local --enable-fortran --enable-fortran2003 --enable-cxx
    $ make
    $ make install

will compile HDF5 and install it to subdirectories in $HOME/local.  By default this will
use the GCC compiler suite; other compilers can be used by setting the CC, CXX and F77
environment variables.  Note the use of ``--enable-fortran2003``; the Fortran 2003
interface is required by HANDE.

pandas
^^^^^^

Pandas can be installed by

.. code-block:: bash

    $ pip install pandas

If you do not have root access, you can install the library locally with:

.. code-block:: bash

    $ pip install pandas --user

Alternatively, where pip is not available, one can install it locally:

.. code-block:: bash

    $ wget https://github.com/pydata/pandas/archive/v0.15.2.tar.gz
    $ tar -xzvf v0.15.2.tar.gz
    $ cd pandas-0.15.2
    $ python setup.py build
    $ python setup.py install

Again, pandas can be installed locally by replacing the final command with:

.. code-block:: bash

    $ python setup.py install --user
