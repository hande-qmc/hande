.. _output_table:

output options
==============

The ``output`` table contains options relating to directing calculation output. This is currently
only compatible with :ref:`fciqmc`, :ref:`ccmc` and :ref:`hilbert`, though extension to other
calculations would be relatively simple.

``filename``
    type: string.

    Optional. Default: 'stdout'.

    Filename to write any calculation output to. If set to default value, all calculation information
    is printed to stdout.
``reprint_sys``
    type: boolean.

    Optional. Deafult: true.

    If true all information on system single particle basis and symmetry that would usually be
    printed during initialisation is also reprinted at the head of any output file. This is useful
    when identifying what system a calculation was performed on a while later.


