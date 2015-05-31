Input file
==========

HANDE is controlled via an input file which is a simple lua script.  This has the
advantage of creating a clean, simple interface to HANDE whilst allowing advanced users to
perform complex simulations without requiring parsing complicated (and perhaps bespoke)
logic in a custom input parser.  Future work will include exposing more of HANDE via
the lua API, thus increasing the flexibility available.

Running a simulation typically involves creating a quantum system (i.e. a collection of
spins/fermions/etc acting under a specified Hamiltonian) and then performing one or more
calculations on that system.  Both tasks involve calling functions from the input file.

The following sections detail options available in each system and calculation
function.  Variables can be required (i.e. must be specified if a function is called) or
optional (in which case the default is stated).  The type (e.g. float, such as 1.234;
integer, such as 1; boolean, either ``true`` or ``false``) of each variable is also given.

Systems
^^^^^^^

All functions which create a system return a pointer to an object (which currently cannot
be manipulated or inspected from lua).  All calculation functions take this system
variable as an argument.

.. toctree::
   :maxdepth: 2

   systems/lattice
   systems/uegs
   systems/generic

Calculations
^^^^^^^^^^^^

.. toctree::
   :maxdepth: 1

   calcs/fci
   calcs/hilbert
   calcs/kinetic_energy
   calcs/fciqmc
   calcs/ccmc
   calcs/dmqmc
   calcs/simple_fciqmc
   calcs/common

Utilities
^^^^^^^^^

.. toctree::
   :maxdepth: 2

   utils/utils

Appendix
^^^^^^^^

.. toctree::
   :maxdepth: 2

   lua_intro
