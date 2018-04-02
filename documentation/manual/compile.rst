.. _compilation:

Compilation
===========

It is possible to configure and build HANDE using CMake or using bare ``make``.

The former requires CMake 3.6 (or newer). It will generate ``Makefile``-s based
on the given configuration parameters and the detected tools and libraries on
your system. It will, in most cases, work out of the box.

The bare ``make`` build offers a higher degree of customisation. Also in this
case a ``Makefile`` will be generated based on a configuration file, of which
you can find examples in the ``config`` folder.

.. toctree::
   :maxdepth: 2
   :glob:

   compile-with-cmake
   compile-with-make
