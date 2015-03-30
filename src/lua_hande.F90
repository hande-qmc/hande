module lua_hande

! A lua interface to HANDE using the AOTUS library.

! A very quick tutorial/summary/notes on interacting with lua:

! 1. Please read the Lua manual on the C API: http://www.lua.org/pil/24.html.
! 2. aotus has many useful wrappers around the Lua C API for interacting with Lua,
!    including transferring information from Fortran to Lua and vice versa.  See
!    lib/aotus/source for the wrapper functions.
! 3. Everything goes via the Lua stack, which is a C pointer (Fortran)/Lua_State* (C).
!    aotus wraps this into a flu_State object, which we should use wherever possible (ie
!    always).
! 4. Note that each function called by lua gets its own private stack.
! 5. The Lua stack is LIFO (last in, first out), so the first argument passed is at the
!    bottom of the stack and the last argument at the top.  Similarly the first argument
!    returned is at the bottom of the stack and the last argument at the top.
! 6. A function can be called from lua if it is 'registered' with the Lua stack and has
!    a very specific interface, e.g.:
!
!        function test_lua_api(L) result(nreturn) bind(c)
!
!            ! Example function callable from a Lua script.
!            ! Must be bind(C) so it can be called from C/Lua.
!            ! Must take a C pointer (which is a Lua stack) as the sole argument.
!            ! Must return an integer which is the number of variables the
!            ! function returns to Lua by pushing to the stack.
!
!            ! In/Out:
!            !    L: lua state (bare C pointer).
!
!            use flu_binding, only: flu_State, flu_copyptr
!            use, intrinsic :: iso_c_binding, only: c_ptr, c_int
!
!            integer(c_int) :: nreturn
!            type(c_ptr), value :: L
!
!            type(flu_State) :: lua_state
!
!            ! Create flu_state from lua state so we can use the nice bindings
!            ! provided by AOTUS.
!            lua_state = flu_copyptr(L)
!
!            ! Number of variables returned on lua stack.
!            nreturn = 0
!
!            ! Get arguments passed to us by lua (if appropriate).
!
!            ! Now do our work...
!            write (6,'(1X,"Hello from fortran!",/)')
!
!        end function test_lua_api
!
! 7. The function can the be called in lua using function_name(arg1, arg2, ...).  In the
!    example above, there are no arguments so it's simply test_lua_api(), assuming
!    test_lua_api is the name provided to lua in the flu_register call.
! 8. Please read the lua documentation and see the examples in the lua_hande* files for
!    more details!

implicit none

contains

    subroutine run_lua_hande(lua_err)

        ! Generic entry point from which a lua script is run.

        ! Out:
        !    lua_err: error code.  Non-zero if there was any problem reading or running
        !          the lua script.

        use aotus_module, only: open_config_file
        use flu_binding, only: flu_State, fluL_newstate, flu_close

        use errors, only: stop_all, warning
        use parallel
        use utils, only: read_file_to_buffer

        integer, intent(out) :: lua_err

        character(255) :: inp_file, err_string
        character(:), allocatable :: buffer
        integer :: buf_len, ierr
        type(flu_State) :: lua_state
        logical :: t_exists

        if (command_argument_count() > 0) then
            ! Input file specified on the command line.
            call get_command_argument(1, inp_file)
            inquire(file=inp_file, exist=t_exists)
            if (.not.t_exists) then
                call stop_all('read_input','File does not exist:'//trim(inp_file))
            end if

            lua_state = fluL_newstate()
            call register_lua_hande_api(lua_state)

            if (parent) then
                write (6,'(a14,/,1X,13("-"),/)') 'Input options'
                call read_file_to_buffer(buffer, inp_file)
                write (6,'(A)') buffer
                buf_len = len(buffer)
                write (6,'(/,1X,13("-"),/)')
            end if

#ifdef PARALLEL
            call mpi_bcast(buf_len, 1, MPI_INTEGER, 0, mpi_comm_world, ierr)
            if (.not.parent) allocate(character(len=buf_len) :: buffer)
            call mpi_bcast(buffer, buf_len, MPI_CHARACTER, 0, mpi_comm_world, ierr)
#endif

            ! Attempt to run script.  If it fails (ie lua_err is non-zero) then try
            ! parsing it as traditional input file for now.
            call open_config_file(lua_state, inp_file, lua_err, err_string)
            ! [todo] - Change to stop once traditional input mode has been removed.
            if (lua_err == 0) then
                call flu_close(lua_state)
            else
                write (6,*) 'aotus/lua error code:', lua_err
                call warning('run_lua_hande', trim(err_string)//'.  Assuming traditional input file...', 2)
            end if

#ifdef PARALLEL
            call mpi_barrier(mpi_comm_world, ierr)
#endif

        else
            ! Maybe read via STDIN?  Only in traditional mode for now...
            ! [todo] - STDIN functionality with lua?
            lua_err = 1
        end if

    end subroutine run_lua_hande

    subroutine register_lua_hande_api(lua_state)

        ! Register HANDE functions which are exposed to Lua scripts with a Lua state
        ! so that they can be called from a Lua script loaded by the same Lua state.

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.

        use flu_binding, only: flu_State, flu_register
        use tests, only: test_lua_api

        type(flu_State), intent(inout) :: lua_state

        call flu_register(lua_state, 'test_lua_api', test_lua_api)

        ! Utilities
        call flu_register(lua_state, 'mpi_root', mpi_root)

        ! Systems
        call flu_register(lua_state, 'hubbard_k', lua_hubbard_k)
        call flu_register(lua_state, 'hubbard_real', lua_hubbard_real)
        call flu_register(lua_state, 'chung_landau', lua_chung_landau)
        call flu_register(lua_state, 'read_in', lua_read_in)
        call flu_register(lua_state, 'heisenberg', lua_heisenberg)
        call flu_register(lua_state, 'ueg', lua_ueg)

        ! Calculations
        call flu_register(lua_state, 'hilbert_space', lua_hilbert_space)
        call flu_register(lua_state, 'kinetic_energy', lua_kinetic_energy)

    end subroutine register_lua_hande_api

    ! --- Helper functions : lua ---

    function mpi_root(L) result(nreturn) bind(c)

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        !    mpi_root()
        ! Returns:
        !    True if on the MPI root processor.

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int
        use parallel, only: parent
        use flu_binding, only: flu_State, flu_copyptr, flu_pushboolean

        integer(c_int) :: nreturn
        type(c_ptr), value :: L
        type(flu_State) :: lua_state

        lua_state = flu_copyptr(L)
        call flu_pushboolean(lua_state, parent)
        nreturn = 1

    end function mpi_root

    ! --- Helper functions : parsing arguments from lua ---

    subroutine warn_unused_args(lua_state, valid_keys, pos)

        ! Print a warning message for the keys not recognised in a table.

        ! In/Out:
        !    lua_state: lua state containing the table.
        ! In:
        !    keys: list of keys recognised in the table.
        !    pos (optional, default 1): position (ie handle) of the table in the
        !        lua stack.

        use flu_binding, only: flu_State, flu_pushnil, flu_next, flu_tolstring, flu_pop, flu_tolstring

        use, intrinsic :: iso_c_binding,  only: c_null_char

        use errors, only: warning

        type(flu_State), intent(inout) :: lua_state
        character(*), intent(in) :: valid_keys(:)
        integer, intent(in), optional :: pos
        integer :: pos_loc, len, j
        character(:), allocatable :: key, key_list
        character, pointer :: str(:)

        pos_loc = 1
        if (present(pos)) pos_loc = pos

        ! Iterate through the table and print key, value pairs.
        ! See example code from the lua api manual: http://pgl.yoyo.org/luai/i/lua_next.
        call flu_pushnil(lua_state)
        do while (flu_next(lua_state, pos_loc))
            ! key is at index -2 and value at index -1
            str => flu_tolstring(lua_state, -2, len)
            allocate(character(size(str)) :: key)
            do j = 1, size(str)
                key(j:j) = str(j)
            end do
            if (all(valid_keys /= key)) then
                if (allocated(key_list)) then
                    key_list = key_list//', '//key
                else
                    key_list = key
                end if
            end if
            ! remove value, keep key for next iteration.
            call flu_pop(lua_state, 1)
            deallocate(key)
        end do

        if (allocated(key_list)) then
            call warning('warn_unused_args', 'The following keywords are not recognised and have been ignored: '//key_list//'.')
        end if

    end subroutine warn_unused_args

    subroutine get_sys_t(lua_state, sys, new)

        ! Get (or create, if necessary) a sys_t object from the lua stack.

        ! This routine assumes that there is a single object on the stack,
        ! which is a table. If this table contains a sys_t object already
        ! (with the key 'sys') then a pointer to this existing object will be
        ! returned. Otherwise a new sys_t object will be created, and a
        ! pointer to this object returned instead.

        ! If it is expected that a sys_t object already exists then the
        ! 'new' variable should not be passed in. A stop_all will then be
        ! called if the sys_t object doesn't exist. If it is not known
        ! if a sys_t object will already exist (or it is expected that it
        ! doesn't) then a 'new' variable should be input.

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        ! Out:
        !    sys: sys_t object from stack/created.
        !    new (optional): If present, set to true if the sys_t object was
        !       allocated rather than passed in.  If not present, an error is
        !       thrown if the 'sys' key is not present.

        use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer
        use errors, only: stop_all
        use flu_binding, only: flu_State, flu_gettop
        use aot_table_module, only: aot_exists, aot_get_val
        use system, only: sys_t

        type(flu_State), intent(inout) :: lua_state
        type(sys_t), pointer, intent(out) :: sys
        logical, optional, intent(out) :: new

        type(c_ptr) :: sys_ptr
        integer :: err

        ! Check if an existing system object was passed in. If so then use
        ! this object instead of creating a new one. thandle=1 because the
        ! table should be the only object on the stack.
        if (aot_exists(lua_state, thandle=1, key='sys')) then
            call aot_get_val(sys_ptr, err, lua_state, thandle=1, key='sys')
            if (err /= 0) call stop_all('get_sys_t', 'Problem receiving sys_t object.')
            call c_f_pointer(sys_ptr, sys)
            if (present(new)) new = .false.
        else
            if (present(new)) then
                new = .true.
            else
                call stop_all('get_sys_t', 'No system object supplied.')
            end if
            allocate(sys)
        end if

    end subroutine get_sys_t

    subroutine set_common_sys_options(lua_state, sys, opts)

        ! Parse system settings common to all (or almost all) system definitions.

        ! In:
        !    opts: handle to the opts table (the main argument passed to the system wrappers).
        ! In/Out:
        !    lua_state: flu/Lua state which has the opts table at the top of the stack.
        !    sys: system (sys_t) object.  On exit the electron number is set, along with symmetry
        !         and spin indices (currently in calc).

        use aot_table_module, only: aot_get_val
        use aot_vector_module, only: aot_get_val
        use errors, only: stop_all
        use flu_binding, only: flu_State
        use system, only: sys_t

        use calc, only: sym_in, ms_in

        type(flu_State), intent(inout) :: lua_state
        type(sys_t), intent(inout) :: sys
        integer, intent(in) :: opts

        integer, allocatable :: cas(:), err_arr(:)
        integer :: err

        call aot_get_val(sys%nel, err, lua_state, opts, 'electrons')
        call aot_get_val(sys%nel, err, lua_state, opts, 'nel')
        call aot_get_val(ms_in, err, lua_state, opts, 'ms')
        call aot_get_val(sym_in, err, lua_state, opts, 'sym')

        call aot_get_val(cas, err_arr, 2, lua_state, opts, key='CAS')
        ! AOTUS returns a vector of size 0 to denote a non-existent vector.
        if (size(cas) == 0) deallocate(cas)
        if (allocated(cas)) then
            if (size(cas) /= 2) call stop_all('set_common_sys_options', &
                            'The CAS option should provide exactly 2 parameters (in an array).')
            sys%CAS = cas
        end if

    end subroutine set_common_sys_options

    subroutine get_ktwist(lua_state, sys, opts)

        ! In:
        !    opts: handle to the opts table (the main argument passed to the system wrappers).
        ! In/Out:
        !    lua_state: flu/Lua state which has the opts table at the top of the stack.
        !    sys: system (sys_t) object.  On entry, sys%lattice%ndim must be set.  On exit the
        !         sys%k_lattice%ktwist array is set, if the option is specified.

        use flu_binding, only: flu_State
        use aot_vector_module, only: aot_get_val

        use const, only: p
        use errors, only: stop_all
        use system, only: sys_t

        type(flu_State), intent(inout) :: lua_state
        type(sys_t), intent(inout) :: sys
        integer, intent(in) :: opts

        real(p), allocatable :: tmp(:)
        integer, allocatable :: err_arr(:)

        call aot_get_val(tmp, err_arr, 3, lua_state, opts, key='twist')
        if (size(tmp) > 0) then
            if (size(tmp) /= sys%lattice%ndim) &
                call stop_all('get_ktwist', 'twist vector not consistent with the dim parameter.')
            allocate(sys%k_lattice%ktwist(sys%lattice%ndim))
            sys%k_lattice%ktwist = tmp
        end if
        deallocate(tmp)

    end subroutine get_ktwist

    subroutine get_lattice(lua_state, sys, opts)

        ! In:
        !    opts: handle to the opts table (the main argument passed to the system wrappers).
        ! In/Out:
        !    lua_state: flu/Lua state which has the opts table at the top of the stack.
        !    sys: system (sys_t) object.  On exit the sys%lattice%ndim and sys%lattice%lattice are set.

        ! NOTE: if get_lattice is called, it is assumed that the lattice option is required.

        use flu_binding, only: flu_State
        use aot_table_ops_module, only: aot_table_open, aot_table_close
        use aot_vector_module, only: aot_get_val

        use const, only: p
        use errors, only: stop_all
        use system, only: sys_t

        type(flu_State), intent(inout) :: lua_state
        type(sys_t), intent(inout) :: sys
        integer, intent(in) :: opts

        real(p), allocatable :: tmp(:)
        integer, allocatable :: err_arr(:)
        integer :: lattice, i

        call aot_table_open(lua_state, opts, lattice, 'lattice')
        if (lattice == 0) call stop_all('get_lattice', 'No lattice specified.')

        call aot_get_val(tmp, err_arr, 3, lua_state, lattice, pos=1)
        if (size(tmp) > 0 .and. size(tmp) <= 3) then
            sys%lattice%ndim = size(tmp)
            allocate(sys%lattice%lattice(sys%lattice%ndim, sys%lattice%ndim))
            sys%lattice%lattice(:,1) = tmp
        else
            call stop_all('get_lattice', 'Unsupported lattice dimension.')
        end if
        deallocate(tmp)

        do i = 2, sys%lattice%ndim
            call aot_get_val(sys%lattice%lattice(:,i), err_arr, lua_state, lattice, pos=i)
        end do

        call aot_table_close(lua_state, lattice)

    end subroutine get_lattice

    subroutine init_generic_system_basis(sys)

        ! A wrapper for initialsing parts of sys_t common to many/all systems.

        ! In/Out:
        !    sys_t: Only registed with init_system and the relevant basis creation before
        !       entry.  On exit, the basis strings and determinants parsing has been
        !       performed.

        use basis_types, only: init_basis_strings, print_basis_metadata
        use determinants, only: init_determinants
        use determinant_enumeration, only: init_determinant_enumeration
        use excitations, only: init_excitations

        use system, only: sys_t, heisenberg

        type(sys_t), intent(inout) :: sys

        ! [todo] - check_input with lua interface?  Should be able to be made simpler now...
        ! call check_input(sys)

        call init_basis_strings(sys%basis)
        call print_basis_metadata(sys%basis, sys%nel, sys%system == heisenberg)
        call init_determinants(sys)
        call init_determinant_enumeration()

        call init_excitations(sys%basis)

    end subroutine init_generic_system_basis

    ! --- System wrappers ---

    function lua_hubbard_k(L) result(nreturn) bind(c)

        ! Create/modify a Hubbard (k-space basis) system.

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        !    hubbard_k{
        !        sys = sys_old -- New sys_t object is created if not passed.
        !        electrons = N,
        !        lattice = { { ... }, { ... }, ... } -- D D-dimensional vectors.
        !        ms = Ms,
        !        sym = sym_index,
        !        U = U,
        !        t = t,
        !        ktwist = {...},    -- D-dimensional vector.
        !    }

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_loc
        use flu_binding, only: flu_State, flu_copyptr, flu_gettop, flu_pushlightuserdata
        use aot_top_module, only: aot_top_get_val
        use aot_table_module, only: aot_table_top, aot_get_val, aot_exists, aot_table_close

        use system, only: sys_t, hub_k, init_system
        use basis, only: init_model_basis_fns
        use momentum_symmetry, only: init_momentum_symmetry

        integer(c_int) :: nreturn
        type(c_ptr), value :: L
        type(flu_State) :: lua_state

        type(sys_t), pointer :: sys
        integer :: opts
        logical :: new, new_basis
        integer :: err
        character(10), parameter :: keys(9) = [character(10) :: 'sys', 'nel', 'electrons', 'lattice', 'U', 't', &
                                                                'ms', 'sym', 'ktwist']

        lua_state = flu_copyptr(L)
        call get_sys_t(lua_state, sys, new)

        ! get a handle to the table...
        opts = aot_table_top(lua_state)

        sys%system = hub_k

        call set_common_sys_options(lua_state, sys, opts)
        call aot_get_val(sys%hubbard%u, err, lua_state, opts, 'U')
        call aot_get_val(sys%hubbard%t, err, lua_state, opts, 't')
        call get_ktwist(lua_state, sys, opts)

        new_basis = aot_exists(lua_state, opts, 'lattice') .or. new

        if (new_basis) then
            call get_lattice(lua_state, sys, opts)
            ! [todo] - deallocate existing basis info and start afresh.
            call init_system(sys)
            call init_model_basis_fns(sys)
            call init_generic_system_basis(sys)
            call init_momentum_symmetry(sys)
        end if

        call warn_unused_args(lua_state, keys, opts)
        call aot_table_close(lua_state, opts)

        call flu_pushlightuserdata(lua_state, c_loc(sys))
        nreturn = 1

    end function lua_hubbard_k

    function lua_hubbard_real(L) result(nreturn) bind(c)

        ! Create/modify a Hubbard (real-space basis) system.

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        !    hubbard_real{
        !        sys = sys_old -- New sys_t object is created if not passed.
        !        electrons = N,
        !        lattice = { { ... }, { ... }, ... } -- D D-dimensional vectors.
        !        ms = Ms,
        !        U = U,
        !        t = t,
        !        finite = true/false,
        !    }

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_loc
        use flu_binding, only: flu_State, flu_copyptr, flu_gettop, flu_pushlightuserdata
        use aot_top_module, only: aot_top_get_val
        use aot_table_module, only: aot_table_top, aot_get_val, aot_exists, aot_table_close

        use system, only: sys_t, hub_real, init_system
        use basis, only: init_model_basis_fns

        integer(c_int) :: nreturn
        type(c_ptr), value :: L
        type(flu_State) :: lua_state

        type(sys_t), pointer :: sys
        integer :: opts
        logical :: new, new_basis
        integer :: err

        character(10), parameter :: keys(8) = [character(10) :: 'sys', 'nel', 'electrons', 'lattice', 'U', 't', &
                                                                'ms', 'finite']

        lua_state = flu_copyptr(L)
        call get_sys_t(lua_state, sys, new)

        ! get a handle to the table...
        opts = aot_table_top(lua_state)

        sys%system = hub_real

        call set_common_sys_options(lua_state, sys, opts)
        call aot_get_val(sys%hubbard%u, err, lua_state, opts, 'U')
        call aot_get_val(sys%hubbard%t, err, lua_state, opts, 't')
        call aot_get_val(sys%real_lattice%finite_cluster, err, lua_state, opts, 'finite')

        new_basis = aot_exists(lua_state, opts, 'lattice') .or. new

        if (new_basis) then
            call get_lattice(lua_state, sys, opts)
            ! [todo] - deallocate existing basis info and start afresh.
            call init_system(sys)
            call init_model_basis_fns(sys)
            call init_generic_system_basis(sys)
        end if

        call warn_unused_args(lua_state, keys, opts)
        call aot_table_close(lua_state, opts)

        call flu_pushlightuserdata(lua_state, c_loc(sys))
        nreturn = 1

    end function lua_hubbard_real

    function lua_chung_landau(L) result(nreturn) bind(c)

        ! Create/modify a Chung--Landau system.

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        !    chung_landau{
        !        sys = sys_old -- New sys_t object is created if not passed.
        !        electrons = N,
        !        lattice = { { ... }, { ... }, ... } -- D D-dimensional vectors.
        !        U = U,
        !        t = t,
        !        finite = true/false,
        !    }

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int
        use flu_binding, only: flu_State, flu_copyptr
        use aot_table_module, only: aot_table_top, aot_exists, aot_table_close

        use errors, only: warning

        integer(c_int) :: nreturn
        type(c_ptr), value :: L
        type(flu_State) :: lua_state
        integer :: opts

        ! The Hamiltonian used by Chung-Landau is essentially the spin-polarised Hubbard model in real
        ! space with a different diagonal matrix element.  Therefore piggy-back on lua_hubbard_real.

        lua_state = flu_copyptr(L)
        opts = aot_table_top(lua_state)
        if (aot_exists(lua_state, opts, 'ms')) then
            call warning('lua_chung_landau', 'The Chung-Landau Hamiltonian is for spinless fermions.  Ignoring the ms keyword.')
        end if

        ! The Hamiltonian used by Chung-Landau is essentially the spin-polarised Hubbard model in real
        ! space with a different diagonal matrix element.  Therefore piggy-back on lua_hubbard_real.
        nreturn = lua_hubbard_real(L)

    end function lua_chung_landau

    function lua_read_in(L) result(nreturn) bind(c)

        ! Create/modify read-in system.

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        !    read_in{
        !        sys = sys_old -- New sys_t object is created if not passed.
        !        electrons = N,
        !        ms = Ms,
        !        int_file = '...',
        !        dipole_int_file = '...'
        !        Lz = true/false
        !        sym = S,
        !        CAS = {cas1, cas2}
        !    }

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_loc
        use flu_binding, only: flu_State, flu_copyptr, flu_gettop, flu_pushlightuserdata
        use aot_top_module, only: aot_top_get_val
        use aot_table_module, only: aot_table_top, aot_get_val, aot_exists, aot_table_close

        use basis, only: init_model_basis_fns
        use momentum_symmetry, only: init_momentum_symmetry
        use read_in_system, only: read_in_integrals
        use system, only: sys_t, read_in, init_system

        integer(c_int) :: nreturn
        type(c_ptr), value :: L
        type(flu_State) :: lua_state

        type(sys_t), pointer :: sys
        type(c_ptr) :: sys_ptr
        integer :: opts
        logical :: new, new_basis
        integer :: err

        character(10), parameter :: keys(8) = [character(10) :: 'sys', 'nel', 'electrons', 'int_file', 'dipole_int_file', 'Lz', &
                                                                'sym', 'ms']

        lua_state = flu_copyptr(L)
        call get_sys_t(lua_state, sys, new)

        ! Get a handle to the table...
        opts = aot_table_top(lua_state)

        sys%system = read_in

        ! Parse table for options...
        call set_common_sys_options(lua_state, sys, opts)

        call aot_get_val(sys%read_in%fcidump, err, lua_state, opts, 'int_file')
        call aot_get_val(sys%read_in%dipole_int_file, err, lua_state, opts, 'dipole_int_file')
        call aot_get_val(sys%read_in%useLz, err, lua_state, opts, 'Lz')

        new_basis = new .or. aot_exists(lua_state, opts, 'int_file') &
                        .or. aot_exists(lua_state, opts, 'CAS')

        if (new_basis) then
            ! [todo] - deallocate existing basis info and start afresh.

            call init_system(sys)
            call read_in_integrals(sys, cas_info=sys%cas)
            call init_generic_system_basis(sys)
        end if

        call warn_unused_args(lua_state, keys, opts)
        call aot_table_close(lua_state, opts)

        sys_ptr = c_loc(sys)
        call flu_pushlightuserdata(lua_state, sys_ptr)
        nreturn = 1

    end function lua_read_in

    function lua_heisenberg(L) result(nreturn) bind(c)

        ! Create/modify a Heisenberg system.

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        !    heisenberg{
        !        sys = sys_old -- New sys_t object is created if not passed.
        !        lattice = { { ... }, { ... }, ... } -- D D-dimensional vectors.
        !        ms = Ms,
        !        J = J,
        !    }

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_loc
        use flu_binding, only: flu_State, flu_copyptr, flu_gettop, flu_pushlightuserdata
        use aot_top_module, only: aot_top_get_val
        use aot_table_module, only: aot_table_top, aot_get_val, aot_exists, aot_table_close

        use system, only: sys_t, heisenberg, init_system
        use basis, only: init_model_basis_fns

        integer(c_int) :: nreturn
        type(c_ptr), value :: L
        type(flu_State) :: lua_state

        type(sys_t), pointer :: sys
        integer :: opts
        logical :: new, new_basis
        integer :: err

        character(10), parameter :: keys(6) = [character(10) :: 'sys', 'ms', 'J', 'lattice', 'magnetic_field', 'staggered_field' ]

        lua_state = flu_copyptr(L)
        call get_sys_t(lua_state, sys, new)

        ! get a handle to the table...
        opts = aot_table_top(lua_state)

        sys%system = heisenberg

        call set_common_sys_options(lua_state, sys, opts)
        call aot_get_val(sys%heisenberg%J, err, lua_state, opts, 'J')
        call aot_get_val(sys%heisenberg%magnetic_field, err, lua_state, opts, 'magnetic_field')
        call aot_get_val(sys%heisenberg%staggered_magnetic_field, err, lua_state, opts, 'staggered_field')

        new_basis = aot_exists(lua_state, opts, 'lattice') .or. new

        if (new_basis) then
            call get_lattice(lua_state, sys, opts)
            ! [todo] - deallocate existing basis info and start afresh.
            call init_system(sys)
            call init_model_basis_fns(sys)
            call init_generic_system_basis(sys)
        end if

        call warn_unused_args(lua_state, keys, opts)
        call aot_table_close(lua_state, opts)

        call flu_pushlightuserdata(lua_state, c_loc(sys))
        nreturn = 1

    end function lua_heisenberg

    function lua_ueg(L) result(nreturn) bind(c)

        ! Create/modify a UEG system.

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        !    ueg{
        !        sys = sys_old -- New sys_t object is created if not passed.
        !        electrons = N,
        !        dim = D,           -- default: 3
        !        rs = density,
        !        cutoff = ecutoff,
        !        ms = Ms,
        !        sym = sym_index,
        !        ktwist = {...},    -- D-dimensional vector
        !        chem_pot = cp
        !    }
        !    Returns: sys_t object.

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_loc
        use flu_binding, only: flu_State, flu_copyptr, flu_gettop, flu_pushlightuserdata
        use aot_top_module, only: aot_top_get_val
        use aot_table_module, only: aot_table_top, aot_get_val, aot_exists, aot_table_close

        use system, only: sys_t, ueg, init_system
        use basis, only: init_model_basis_fns
        use momentum_symmetry, only: init_momentum_symmetry
        use ueg_system, only: init_ueg_proc_pointers

        integer(c_int) :: nreturn
        type(c_ptr), value :: L
        type(flu_State) :: lua_state

        type(sys_t), pointer :: sys
        type(c_ptr) :: sys_ptr
        integer :: opts
        logical :: new_basis, new
        integer :: err

        character(10), parameter :: keys(10) = [character(10) :: 'sys', 'cutoff', 'dim', 'rs', 'nel', 'electrons', &
                                               'ms', 'sym', 'ktwist', 'chem_pot']

        lua_state = flu_copyptr(L)
        call get_sys_t(lua_state, sys, new)

        ! Get a handle to the table...
        opts = aot_table_top(lua_state)

        sys%system = ueg
        sys%lattice%ndim = 3

        ! Parse table for options...
        call set_common_sys_options(lua_state, sys, opts)

        new_basis = aot_exists(lua_state, opts, 'cutoff') .or. &
                    aot_exists(lua_state, opts, 'dim')    .or. &
                    aot_exists(lua_state, opts, 'rs')    .or. &
                    aot_exists(lua_state, opts, 'nel')    .or. &
                    aot_exists(lua_state, opts, 'electrons') .or. new

        call aot_get_val(sys%ueg%ecutoff, err, lua_state, opts, 'cutoff')
        call aot_get_val(sys%ueg%r_s, err, lua_state, opts, 'rs')
        call aot_get_val(sys%ueg%chem_pot, err, lua_state, opts, 'chem_pot')
        call aot_get_val(sys%lattice%ndim, err, lua_state, opts, 'dim')

        call get_ktwist(lua_state, sys, opts)

        if (new_basis) then
            ! [todo] - deallocate existing basis info and start afresh.

            call init_system(sys)
            call init_model_basis_fns(sys)
            call init_generic_system_basis(sys)
            call init_momentum_symmetry(sys)
            call init_ueg_proc_pointers(sys%lattice%ndim, sys%ueg)
        end if

        call warn_unused_args(lua_state, keys, opts)
        call aot_table_close(lua_state, opts)

        sys_ptr = c_loc(sys)
        call flu_pushlightuserdata(lua_state, sys_ptr)
        nreturn = 1

    end function lua_ueg

    ! --- Calculation wrappers ---

    function lua_hilbert_space(L) result(nresult) bind(c)

        ! Run a Monte Carlo calculation to estimate the size of the Hilbert
        ! space.

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        !    hilbert_space{
        !           sys = sys_t,
        !           ncycles = n,
        !           ex_level = level,
        !           rng_seed = seed,
        !           reference = { ... },  -- nel-dimensional vector.
        !    }

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_f_pointer

        use flu_binding, only: flu_State, flu_copyptr
        use aot_top_module, only: aot_top_get_val
        use aot_table_module, only: aot_table_top, aot_get_val, aot_exists, aot_table_close
        use aot_vector_module, only: aot_get_val

        use errors, only: stop_all
        use system, only: sys_t

        use hilbert_space, only: estimate_hilbert_space

        integer(c_int) :: nresult
        type(c_ptr), value :: L

        type(flu_State) :: lua_state

        type(c_ptr) :: sys_ptr
        type(sys_t), pointer :: sys
        integer :: truncation_level, ncycles, rng_seed
        integer, allocatable :: ref_det(:)
        integer :: opts, err
        integer, allocatable :: err_arr(:)
        logical :: have_seed
        character(10), parameter :: keys(5) = [character(10) :: 'sys', 'ncycles', 'ex_level', 'reference', 'rng_seed']

        lua_state = flu_copyptr(L)
        call get_sys_t(lua_state, sys)

        opts = aot_table_top(lua_state)

        call aot_get_val(ncycles, err, lua_state, opts, 'ncycles')
        if (err /= 0) call stop_all('lua_hilbert_space', 'Number of cycles not supplied.')
        call aot_get_val(truncation_level, err, lua_state, opts, 'ex_level', default=-1)
        call aot_get_val(ref_det, err_arr, sys%nel, lua_state, opts, key='reference')
        have_seed = aot_exists(lua_state, opts, 'rng_seed')
        call aot_get_val(rng_seed, err, lua_state, opts, 'rng_seed')

        call warn_unused_args(lua_state, keys, opts)
        call aot_table_close(lua_state, opts)

        ! AOTUS returns a vector of size 0 to denote a non-existent vector.
        if (size(ref_det) == 0) deallocate(ref_det)
        if (allocated(ref_det)) then
            if (size(ref_det) /= sys%nel) call stop_all('lua_hilbert_space', &
                            'Reference determinant does not match the number of electrons in system.')
        end if

        if (have_seed) then
            call estimate_hilbert_space(sys, truncation_level, ncycles, ref_det, rng_seed)
        else
            call estimate_hilbert_space(sys, truncation_level, ncycles, ref_det)
        end if

        ! [todo] - return estimate of space and error to lua.
        nresult = 0

    end function lua_hilbert_space

    function lua_kinetic_energy(L) result(nresult) bind(c)

        ! Run a Monte Carlo calculation to estimate the kinetic energy of a
        ! system in the canonical ensemble.

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        !    kinetic_energy{
        !           sys = sys_t,
        !           nattempts = npop,
        !           ncycles = n,
        !           beta = beta,
        !           fermi_temperature = true/false,
        !           rng_seed = seed,
        !    }

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_f_pointer
        use const, only: p

        use flu_binding, only: flu_State, flu_copyptr
        use aot_top_module, only: aot_top_get_val
        use aot_table_module, only: aot_table_top, aot_get_val, aot_exists, aot_table_close
        use aot_vector_module, only: aot_get_val

        use calc, only: ms_in
        use errors, only: stop_all
        use system, only: sys_t

        use canonical_kinetic_energy, only: estimate_kinetic_energy

        integer(c_int) :: nresult
        type(c_ptr), value :: L

        type(flu_State) :: lua_state

        type(c_ptr) :: sys_ptr
        type(sys_t), pointer :: sys
        integer :: opts, err, seed, ncycles, nattempts
        real(p) :: beta
        logical :: fermi_temperature, have_seed
        character(10), parameter :: keys(6) = [character(10) :: 'sys', 'nattempts', 'ncycles', 'beta', &
                                                                'fermi_temperature', 'rng_seed']

        lua_state = flu_copyptr(L)
        call get_sys_t(lua_state, sys)

        opts = aot_table_top(lua_state)

        call aot_get_val(nattempts, err, lua_state, opts, 'nattempts')
        if (err /= 0) call stop_all('lua_kinetic_energy', 'nattempts: Number of attempts/cycle not supplied.')
        call aot_get_val(ncycles, err, lua_state, opts, 'ncycles')
        if (err /= 0) call stop_all('lua_kinetic_energy', 'ncycles: Number of cycles not supplied.')

        call aot_get_val(beta, err, lua_state, opts, 'beta')
        if (err /= 0) call stop_all('lua_kinetic_energy', 'beta: target temperature not supplied.')
        call aot_get_val(fermi_temperature, err, lua_state, opts, 'fermi_temperature', default=.false.)

        if (sys%ueg%chem_pot == huge(1.0_p)) call stop_all('lua_kinetic_energy', &
                                                           'chem_pot: chemical potential not supplied.')

        have_seed = aot_exists(lua_state, opts, 'rng_seed')
        call aot_get_val(seed, err, lua_state, opts, 'rng_seed')

        call warn_unused_args(lua_state, keys, opts)
        call aot_table_close(lua_state, opts)

        if (have_seed) then
            call estimate_kinetic_energy(sys, fermi_temperature, beta, nattempts, ncycles, seed)
        else
            call estimate_kinetic_energy(sys, fermi_temperature, beta, nattempts, ncycles)
        end if

        ! [todo] - return estimate of kinetic energy and error to lua.
        nresult = 0

    end function lua_kinetic_energy

end module lua_hande
