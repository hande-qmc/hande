.. _fciqmc_tutorial:

Full Configuration Interaction Quantum Monte Carlo
==================================================

In this tutorial we will run FCIQMC on the 18-site 2D Hubbard model at half filling.

First, we will set up the system and estimate the number of determinants in Hilbert space
with the desired symmetry using a Monte Carlo approach.

We are interested in the state with zero crystal momentum, as there is theoretical work
showing this will be the symmetry of the overall ground state.  HANDE uses an indexing
scheme for the symmetry label.  The easiest way to find this out is to run an input file
which only contains the system definition: 

.. code-block:: lua

    hubbard = hubbard_k {
        lattice = {
                    { 3,  3 },
                    { 3, -3 },
                  },
        electrons = 18, 
        ms = 0,
    }

Running this using

.. code-block:: bash

    $ hande.x hubbard.lua > hubbard.out

The output contains a symmetry table which informs us that the wavevector :math:`(0,0)`
corresponds to the index 1.
        

.. Create system
.. Run calculation -> find plateau
.. Run production calculation
.. Analyse data
