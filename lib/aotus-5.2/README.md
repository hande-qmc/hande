NOTE::

We only include files from AOTUS needed for HANDE.  Please see upstream for the complete
distribution.  In particular, tests, examples, documentation, waf build files and extended double and
quadruple precision source files are not included.

Advanced Options and Tables in Universal Scripting
==================================================

The AOTUS library provides a Fortran wrapper around the C-API of the
[Lua](http://www.lua.org) scripting language, allowing a convenient usage of Lua
scripts as configuration files in Fortran applications.
Please have a look at the [Wiki](https://bitbucket.org/haraldkl/aotus/wiki/Home)
for more information on its usage.
And for a detailed interface reference visit its
[Doxygen generated documentation](https://geb.sts.nt.uni-siegen.de/doxy/aotus)

*This library is released under a simplified MIT licence, please have a look into the COPYRIGHT file for details.*

Aotus is part of the APES suite, for which there is a
[mailing list](https://listserv.uni-siegen.de/cgi-bin/mailman/listinfo/apes)
where questions can be asked.


How To Build
------------

[Waf](http://code.google.com/p/waf/) is used as build system.
Run:

~~~~~~~~~~~{.sh}
./waf configure build
~~~~~~~~~~~

to build the aotus library.
If you want to select a specific Fortran compiler, set the environment variable
*FC*.
And for a specific C compiler, set the environment variable *CC*.
The Fortran compiler flags are set with the help of fc_flags, which provide
a set of compiler flag combinations for various compilers.
They are found in the fc_flags.py file in the root directory of the project.

By running:

~~~~~~~~~~~{.sh}
./waf --help
~~~~~~~~~~~

you get a list of available options to the waf script.


What is Built
-------------

For your convenience the Lua library is included in version 5.2.3 (released
2013-12-07).
Its objects are completely gathered into the final *libaotus* library, so it is
only necessary to link against this single static library to gain the
configuration features of aotus in your Fortran application.
Due to the compiler specific module information required by any application
using the libaotus, the suggested approach to incorporate libaotus is to include
its building in the build process of the final application. This is straight
forward if waf is used for the complete project. But also in other build
environments it should not be too hard to make use of the generated *build*
directory.
Yet if you would rather install the *libaotus.a* and the module files into a
*$PREFIX* directory, you can make use of:

~~~~~~~~~~~{.sh}
./waf install
~~~~~~~~~~~

The default build process will also create some unit test executables and
execute them to ensure functionality of the various parts in the library.

The doxygen documentation can be built by running:

~~~~~~~~~~~{.sh}
./waf doxy
~~~~~~~~~~~

This will build a html directory in the build directory with the resulting
documentation. Note, that this requires an installed doxygen.
It is also online available at
[Aotus documentation](https://geb.sts.nt.uni-siegen.de/doxy/aotus).

### Example

There is an example program built, called aotus_sample, which you will find in
the *build* directory.
It can be used with the provided *config.lua* in the *sample* directory, where
also the source of this small program is found.

Related Projects
----------------

Some projects with similar goals or related information:

* [f2k3-lua](https://github.com/MaikBeckmann/f2k3-lua/tree/simple)
* [FortLua](https://github.com/adolgert/FortLua)
