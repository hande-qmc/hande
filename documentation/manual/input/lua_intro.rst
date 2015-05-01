A short introduction to lua
===========================

Lua is a lightweight programming language which is easy to embed and is well-suited to the
task of controlling a simulation.  For a quick introduction to lua, please read
`Learn Lua in 15 Minutes <http://tylerneylon.com/a/learn-lua/>`_.  However, for most cases
the input file format can be treated as follows:

Assignment is perfomed by setting a variable name equal to an object, e.g.

.. code-block:: lua

    pi = 3.141592654

Strings are created by enclosing characters in quotation marks:

.. code-block:: lua

    msg = 'hello world'
    
and boolean variables can be set using the ``true`` and ``false`` keywords:

.. code-block:: lua

    yes = true
    yes = false

A key data structure in lua is the `table`, which serves both as an array and an
associative array or map, and is denoted using braces.  First, the following creates a table
to hold a 1D vector:

.. code-block:: lua

    v = { 1, 2, 3 }

whilst using key=value pairs creates a table as an associative array:

.. code-block:: lua

    v = { x = 3, y = 4, type = 'dual' }

Tables can be nested.

Functions are called using:

.. code-block:: lua

    x = fname(arg1, arg2, ...)

where `fname` is the name of the function, which returns a single value (which is
stored in `x` in the above example).  Keywords can be passed in by using a table.  If the
function takes a single table, then the parentheses need not be included, such that the
following calls are identical:

.. code-block:: lua

    x = fname1({ x = 3, y = 4, type = 'dual'})
    x = fname1{ x = 3, y = 4, type = 'dual'}

All options are passed into HANDE by using a table as an associative array.  Each function
exposed by HANDE to the lua script takes a single (nested) table.

.. warning::

    lua, and by extension the HANDE input file, is case sensitive.
