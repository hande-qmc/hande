Debugging options
-----------------

There are a couple of compilation options to help with debugging HANDE.

* The ``-g`` option to ``tools/mkconfig.py`` enables compiler options
  for warnings and run-time checking.

* The ``-DDEBUG`` preprocessor flag enables additional debugging output.
  Currently this is stack traces when ``stop_all`` is called to terminate
  with an error - the addresses given can be converted to ``file:line number``
  information with ``addr2line``. For example:

  .. code-block:: bash

      $ /path/to/hande.x test.lua
      [...]
      /usr/lib/libasan.so.4(+0x55c60)[0x7fc1c5000c60]
      ./bin/hande.x(+0x609e80)[0x5614f6742e80]
      ./bin/hande.x(+0x469177)[0x5614f65a2177]
      ./bin/hande.x(+0x1aa839)[0x5614f62e3839]
      ./bin/hande.x(+0x1aa8f1)[0x5614f62e38f1]
      /usr/lib/libc.so.6(__libc_start_main+0xea)[0x7fc1c17e1f4a]
      ./bin/hande.x(+0xa93ca)[0x5614f61e23ca]

      ERROR.
      HANDE stops in subroutine: run_hande_lua.
      Reason: File does not exist:test.lua
      EXITING...

  The addresses are not conserved between builds, compilers, optimisation levels and so
  on. The filenames can be included if ``-rdynamic`` is included in the linker flags::

       /usr/lib/libasan.so.4(+0x55c60)[0x7f1943b93c60]
       ./bin/hande.x(__errors_MOD_stop_all+0x1c2)[0x55bf6aae5d30]
       ./bin/hande.x(__lua_hande_MOD_run_lua_hande+0x6e1)[0x55bf6a945027]
       ./bin/hande.x(+0x1c66e9)[0x55bf6a6866e9]
       ./bin/hande.x(main+0x36)[0x55bf6a6867a1]
       /usr/lib/libc.so.6(__libc_start_main+0xea)[0x7f1940374f4a]
       ./bin/hande.x(_start+0x2a)[0x55bf6a58527a]
       
       ERROR.
       HANDE stops in subroutine: run_hande_lua.
       Reason: File does not exist:test.lua
       EXITING...

  The actual source (file and linenumber) can be found using either ``addr2line`` or
  ``gdb`` (easier). With addr2line:

  .. code-block:: bash

      $ addr2line +0x469177 -e /path/to/hande.x
      /home/james/hande/src/hande-bug-fix/src/lua_hande.F90:92
  
  and with gdb:

  .. code-block:: bash

      $ gdb /path/to/hande.x
      (gdb) list *__lua_hande_MOD_run_lua_hande+0x6e1
      0x485027 is in lua_hande::run_lua_hande (src/lua_hande.F90:92).
      87	
      88	            ! Read input file on parent and broadcast to all other processors.
      89	            if (parent) then
      90	                call get_command_argument(1, inp_file)
      91	                inquire(file=inp_file, exist=t_exists)
      92	                if (.not.t_exists) call stop_all('run_hande_lua','File does not exist:'//trim(inp_file))
      93	
      94	                write (6,'(a14,/,1X,13("-"),/)') 'Input options'
      95	                call read_file_to_buffer(buffer, inp_file)
      96	                write (6,'(A)') trim(buffer)

  where ``__lua_hande_MOD_run_lua_hande+0x6e1`` was the address of interest from the
  stacktrace. Note the ``*`` prefix. The same can be done with just the bare address if
  ``-rdynamic`` isn't used.

  addr2line is a little more involved if ``-rdynamic`` is used -- one needs to find the
  offset for the function of interest and hence find the line address by adding the
  address relative to the function start. For example, in::

      ./bin/hande.x(__errors_MOD_stop_all+0x1c2)[0x55bf6aae5d30]

  the function is ``__errors_MOD_stop_all`` (``stop_all`` in the ``errors`` module). The
  start address can be found from ``objdump``:

  .. code-block:: bash

      $ objdump -T /path/to/hande.x | grep __errors_MOD_stop_all
      0000000000625b6e g    DF .text	0000000000000c0c  Base        __errors_MOD_stop_all

  The first field is the address in **hexadecimal**. Hence:

  .. code-block:: bash

      $ python -c 'print(hex(0x0000000000625b6e+0x1c2))'
      0x625d30
      $ addr2line -e /path/to/hande.x 0x625d30
      /path/to/hande/lib/local/error_handling.F90:60

  Note the explicit ``0x`` prefix for the start address of ``__errors_MOD_stop_all``.
