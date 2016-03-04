Debugging options
-----------------

There are a couple of compilation options to help with debugging HANDE.

* The ``-g`` option to ``tools/mkconfig.py`` enables compiler options
  for warnings and run-time checking.

* The ``-DDEBUG`` preprocessor flag enables additional debugging output.
  Currently this is stack traces when ``stop_all`` is called to terminate
  with an error - the addresses given can be converted to ``file:line number``
  information with ``addr2line``.
