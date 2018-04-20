module lua_hande_system

! Lua wrappers to create/modify systems.  See top-level comments in lua_hande
! for more details about working with the Lua API.

implicit none

contains

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
        use parallel, only: parent
        use errors, only: stop_all
        use lua_hande_utils, only: warn_unused_args
        use flu_binding, only: flu_State, flu_gettop
        use aot_table_module, only: aot_exists, aot_get_val, aot_table_top
        use aot_table_ops_module, only: aot_table_open, aot_table_close
        use system, only: sys_t
        use lua_hande_utils, only: get_userdata

        type(flu_State), intent(inout) :: lua_state
        type(sys_t), pointer, intent(out) :: sys
        logical, optional, intent(out) :: new

        type(c_ptr) :: sys_ptr
        integer :: table, sys_table
        logical :: have_sys_entry

        table = aot_table_top(lua_state)
        if (table /= 0) then
            have_sys_entry = .false.
            have_sys_entry = aot_exists(lua_state, thandle=table, key='sys')
        else
            have_sys_entry = .false.
        end if

        ! Check if an existing system object was passed in. If so then use
        ! this object instead of creating a new one.
        if (have_sys_entry) then
            call aot_table_open(lua_state, table, sys_table, key='sys')
            call get_userdata(lua_state, sys_table, "sys", sys_ptr)
            call aot_table_close(lua_state, sys_table)
            call c_f_pointer(sys_ptr, sys)
            if (present(new)) new = .false.
        else
            if (present(new)) then
                new = .true.
            else if (parent) then
                call stop_all('get_sys_t', 'No system object supplied.')
            end if
            allocate(sys)
        end if

    end subroutine get_sys_t

    subroutine push_sys(lua_state, sys, new)

        ! Add a table containing the passed system to the lua stack for returning

        ! In/Out:
        !   lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !   sys: the sys object to return to lua
        !   new: if false, do not create a new system object

        use, intrinsic :: iso_c_binding, only: c_loc
        use flu_binding, only: flu_State, flu_pushlightuserdata, flu_pushstring, flu_settable, flu_pushcclosure, fluL_setmetatable
        use aot_table_ops_module, only: aot_table_open, aot_table_close

        use system, only: sys_t

        type(sys_t), pointer, intent(in) :: sys
        type(flu_state), intent(inout) :: lua_state
        logical, intent(in) :: new

        integer :: table

        if (new) then
            ! Create table to become sys object
            call aot_table_open(lua_state, thandle=table)

            ! Add sys pointer as t.sys
            call flu_pushstring(lua_state, "sys")
            call flu_pushlightuserdata(lua_state, c_loc(sys))
            call flu_settable(lua_state, table)

            ! Add deallocation function as t:free()
            call flu_pushstring(lua_state, "free")
            call flu_pushcclosure(lua_state, lua_dealloc_sys, 0)
            call flu_settable(lua_state, table)

            ! Set metatable to mark for finalisation.  Note metatable is created in register_lua_hande_api.
            call fluL_setmetatable(lua_state, "sys")
        else
            ! Push already existing sys table/object onto stack
            call aot_table_open(lua_state, 1, table, 'sys')
        end if


    end subroutine push_sys

    subroutine set_common_sys_options(lua_state, sys, verbose, opts)

        ! Parse system settings common to all (or almost all) system definitions.

        ! In:
        !    opts: handle to the opts table (the main argument passed to the system wrappers).
        ! In/Out:
        !    lua_state: flu/Lua state which has the opts table at the top of the stack.
        !    sys: system (sys_t) object.  On exit the electron number is set, along with symmetry
        !         and spin indices (currently in calc).
        ! Out:
        !    verbose: set verbosity of output for basis information.

        use aot_table_module, only: aot_get_val, aot_exists
        use aot_vector_module, only: aot_get_val
        use parallel, only: parent
        use errors, only: stop_all
        use lua_hande_utils, only: warn_unused_args
        use flu_binding, only: flu_State
        use system, only: sys_t

        type(flu_State), intent(inout) :: lua_state
        type(sys_t), intent(inout) :: sys
        logical, intent(out) :: verbose
        integer, intent(in) :: opts

        integer, allocatable :: cas(:), err_arr(:)
        integer :: err
        character(len=10) :: str

        call aot_get_val(sys%nel, err, lua_state, opts, 'electrons')
        call aot_get_val(sys%nel, err, lua_state, opts, 'nel')
        call aot_get_val(sys%Ms, err, lua_state, opts, 'ms')

        call aot_get_val(verbose, err, lua_state, opts, 'verbose', default=.true.)
        ! Treat sym more carefully for various options.
        if (aot_exists(lua_state, opts, 'sym')) then
            call aot_get_val(sys%symmetry, err, lua_state, opts, 'sym')
            if (err /= 0) then
                sys%symmetry = huge(0)
                ! Might have passed string in?
                call aot_get_val(str, err, lua_state, opts, 'sym')
                sys%tot_sym = trim(str) == 'tot_sym'
                sys%aufbau_sym = trim(str) == 'aufbau'
            else
                ! sys%tot_sym is already initialized to .false.
                sys%aufbau_sym = .false.
            end if
            if (.not. (sys%tot_sym .or. sys%aufbau_sym .or. sys%symmetry < huge(0))) &
                call stop_all('set_common_sys_options', 'sym input not recognised. Please &
                    &refer to the documentation for guidelines to valid input.')
        end if

        call aot_get_val(cas, err_arr, 2, lua_state, opts, key='CAS')
        ! AOTUS returns a vector of size 0 to denote a non-existent vector.
        if (size(cas) == 0) deallocate(cas)
        if (allocated(cas)) then
            if (size(cas) /= 2 .and. parent) call stop_all('set_common_sys_options', &
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
        use lua_hande_utils, only: warn_unused_args
        use system, only: sys_t

        type(flu_State), intent(inout) :: lua_state
        type(sys_t), intent(inout) :: sys
        integer, intent(in) :: opts

        real(p), allocatable :: tmp(:)
        integer, allocatable :: err_arr(:)

        call aot_get_val(tmp, err_arr, 3, lua_state, opts, key='twist')
        if (size(tmp) > 0) then
            allocate(sys%k_lattice%ktwist(size(tmp)))
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

        use parallel, only: parent
        use errors, only: stop_all
        use lua_hande_utils, only: warn_unused_args
        use system, only: sys_t

        type(flu_State), intent(inout) :: lua_state
        type(sys_t), intent(inout) :: sys
        integer, intent(in) :: opts

        integer, allocatable :: tmp(:)
        integer, allocatable :: err_arr(:)
        integer :: lattice, i

        call aot_table_open(lua_state, opts, lattice, 'lattice')
        if (lattice == 0 .and. parent) call stop_all('get_lattice', 'No lattice specified.')

        call aot_get_val(tmp, err_arr, 3, lua_state, lattice, pos=1)
        if (size(tmp) > 0 .and. size(tmp) <= 3) then
            sys%lattice%ndim = size(tmp)
            allocate(sys%lattice%lattice(sys%lattice%ndim, sys%lattice%ndim))
            sys%lattice%lattice(:,1) = tmp
        else
            if (parent) call stop_all('get_lattice', 'Unsupported lattice dimension.')
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
        use excitations, only: init_excitations

        use system, only: sys_t, heisenberg

        type(sys_t), intent(inout) :: sys

        call init_basis_strings(sys%basis)
        call print_basis_metadata(sys%basis, sys%nel, sys%system == heisenberg)
        ! NOTE: RAS is currently disabled.
        call init_determinants(sys, sys%nel)

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
        !        twist = {...},    -- D-dimensional vector.
        !    }

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int
        use flu_binding, only: flu_State, flu_copyptr, flu_gettop
        use aot_top_module, only: aot_top_get_val
        use aot_table_module, only: aot_table_top, aot_get_val, aot_exists, aot_table_close

        use lua_hande_utils, only: warn_unused_args, register_timing
        use system, only: sys_t, hub_k, init_system
        use calc_system_init, only: set_spin_polarisation
        use basis, only: init_model_basis_fns
        use momentum_symmetry, only: init_momentum_symmetry
        use check_input, only: check_sys
        use parallel, only: parent
        use errors, only: stop_all

        integer(c_int) :: nreturn
        type(c_ptr), value :: L
        type(flu_State) :: lua_state

        type(sys_t), pointer :: sys
        integer :: opts
        real :: t1, t2
        logical :: new, new_basis, verbose
        integer :: err
        character(10), parameter :: keys(10) = [character(10) :: 'sys', 'nel', 'electrons', 'lattice', 'U', 't', &
                                                                'ms', 'sym', 'twist', 'verbose']

        call cpu_time(t1)

        lua_state = flu_copyptr(L)
        call get_sys_t(lua_state, sys, new)

        ! get a handle to the table...
        opts = aot_table_top(lua_state)

        sys%system = hub_k

        call set_common_sys_options(lua_state, sys, verbose, opts)
        call aot_get_val(sys%hubbard%u, err, lua_state, opts, 'U')
        call aot_get_val(sys%hubbard%t, err, lua_state, opts, 't')
        call get_ktwist(lua_state, sys, opts)

        new_basis = aot_exists(lua_state, opts, 'lattice') .or. new

        if (new_basis) then
            call get_lattice(lua_state, sys, opts)
            ! [todo] - deallocate existing basis info and start afresh.
            call init_system(sys)
            if (parent) call check_sys(sys)
            call init_model_basis_fns(sys, verbose=verbose)
            call init_generic_system_basis(sys)
            call init_momentum_symmetry(sys)
        end if

        if (parent) then
            ! This check can't go with the rest of check_sys as it has to happen after the basis is initialised.
            if (sys%symmetry < huge(0) .and. (sys%symmetry > sys%sym_max .or. sys%symmetry < sys%sym0)) &
                call stop_all('lua_hubbard_k', 'Symmetry out of range.')
        end if

        call warn_unused_args(lua_state, keys, opts)

        call push_sys(lua_state, sys, new_basis)
        nreturn = 1

        call cpu_time(t2)
        call register_timing(lua_state, "Momentum space hubbard model initialisation", t2-t1)

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

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int
        use system, only: hub_real

        integer(c_int) :: nreturn
        type(c_ptr), value :: L

        nreturn = real_lattice_wrapper(L, hub_real)

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
        use system, only: chung_landau

        integer(c_int) :: nreturn
        type(c_ptr), value :: L

        ! The Hamiltonian used by Chung-Landau is essentially the spin-polarised Hubbard model in real
        ! space with a different diagonal matrix element.  Therefore piggy-back on lua_hubbard_real.
        nreturn = real_lattice_wrapper(L, chung_landau)

    end function lua_chung_landau

    function real_lattice_wrapper(L, system_type) result(nreturn)

        ! Create/modify a Chung--Landau or Hubbard (real-space basis) system.

        ! In/Out:
        !    L: lua state (bare C pointer).
        ! In:
        !    system_type: system type to create/modify.

        ! See comments for lua_hubbard_real and lua_chung_landau for accepted keys.

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int
        use flu_binding, only: flu_State, flu_copyptr, flu_gettop
        use aot_top_module, only: aot_top_get_val
        use aot_table_module, only: aot_table_top, aot_get_val, aot_exists, aot_table_close

        use lua_hande_utils, only: warn_unused_args, register_timing
        use system, only: sys_t, chung_landau, hub_real, init_system
        use basis, only: init_model_basis_fns
        use real_lattice, only: init_real_space
        use check_input, only: check_sys
        use parallel, only: parent
        use errors, only: stop_all

        integer(c_int) :: nreturn
        type(c_ptr) :: L
        integer, intent(in) :: system_type
        type(flu_State) :: lua_state

        type(sys_t), pointer :: sys
        integer :: opts
        real :: t1, t2
        logical :: new, new_basis, verbose
        integer :: err

        character(10), parameter :: keys(9) = [character(10) :: 'sys', 'nel', 'electrons', 'lattice', 'U', 't', &
                                                                'finite', 'ms', 'verbose']

        call cpu_time(t1)

        lua_state = flu_copyptr(L)
        call get_sys_t(lua_state, sys, new)

        ! get a handle to the table...
        opts = aot_table_top(lua_state)

        sys%system = system_type

        call set_common_sys_options(lua_state, sys, verbose, opts)
        call aot_get_val(sys%hubbard%u, err, lua_state, opts, 'U')
        call aot_get_val(sys%hubbard%t, err, lua_state, opts, 't')
        call aot_get_val(sys%real_lattice%finite_cluster, err, lua_state, opts, 'finite')

        ! Only one symmetry sector, so we can set it here (needed for dmqmc initialisation)
        sys%symmetry = 1

        new_basis = aot_exists(lua_state, opts, 'lattice') .or. new

        if (new_basis) then
            call get_lattice(lua_state, sys, opts)
            ! [todo] - deallocate existing basis info and start afresh.
            call init_system(sys)
            if (parent) call check_sys(sys)
            call init_model_basis_fns(sys, verbose=verbose)
            call init_generic_system_basis(sys)
            call init_real_space(sys)
        end if

        call warn_unused_args(lua_state, keys, opts)

        call push_sys(lua_state, sys, new_basis)
        nreturn = 1

        call cpu_time(t2)
        select case (system_type)
        case (hub_real)
            call register_timing(lua_state, "Real space hubbard model initialisation", t2-t1)
        case (chung_landau)
            call register_timing(lua_state, "Chung-landau initialisation", t2-t1)
        case default
            call stop_all("real_lattice_wrapper", "Unknown system type")
        end select

    end function real_lattice_wrapper

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
        !        complex = true/false,
        !        max_broadcast_chunk = block_size,
        !    }

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int
        use flu_binding, only: flu_State, flu_copyptr, flu_gettop
        use aot_top_module, only: aot_top_get_val
        use aot_table_module, only: aot_table_top, aot_get_val, aot_exists, aot_table_close

        use lua_hande_utils, only: warn_unused_args, register_timing
        use hdf5_helper, only: ishdf5_wrapper
        use point_group_symmetry, only: print_pg_symmetry_info
        use momentum_sym_read_in, only: print_mom_sym_info
        use read_in_system, only: read_in_integrals
        use system, only: sys_t, read_in, init_system
        use check_input, only: check_sys
        use parallel, only: parent
        use errors, only: stop_all
        use hdf5_system, only: read_system_hdf5

        integer(c_int) :: nreturn
        type(c_ptr), value :: L
        type(flu_State) :: lua_state

        type(sys_t), pointer :: sys
        integer :: opts
        real :: t1, t2
        logical :: new, new_basis, verbose, hdf5, t_exists
        integer :: err

        character(20), parameter :: keys(12) = [character(20) :: 'sys', 'nel', 'electrons', 'int_file', 'dipole_int_file', 'Lz', &
                                                                'sym', 'ms', 'CAS', 'complex', 'verbose', 'max_broadcast_chunk']

        call cpu_time(t1)

        lua_state = flu_copyptr(L)
        call get_sys_t(lua_state, sys, new)

        ! Get a handle to the table...
        opts = aot_table_top(lua_state)

        sys%system = read_in

        ! Parse table for options...
        call set_common_sys_options(lua_state, sys, verbose, opts)

        call aot_get_val(sys%read_in%fcidump, err, lua_state, opts, 'int_file')

        call aot_get_val(sys%read_in%dipole_int_file, err, lua_state, opts, 'dipole_int_file')

        call aot_get_val(sys%read_in%max_broadcast_chunk, err, lua_state, opts, 'max_broadcast_chunk')


        if (parent) then
            ! Verify that the specified file exists and check whether it is HDF5 or text.
            inquire(file=sys%read_in%fcidump, exist=t_exists)
            if (.not.t_exists) call stop_all('lua_read_in', 'FCIDUMP does not exist: '//trim(sys%read_in%fcidump))
        end if
        call ishdf5_wrapper(sys%read_in%fcidump, hdf5, err)

        if (hdf5) then
            call init_system(sys)
            call read_system_hdf5(sys, verbose)
            if (parent) call check_sys(sys)
            new_basis = .true.
            if (sys%momentum_space) then
                call print_mom_sym_info(sys)
            else
                call print_pg_symmetry_info(sys)
            end if
        else

            call aot_get_val(sys%read_in%useLz, err, lua_state, opts, 'Lz')
            call aot_get_val(sys%read_in%comp, err, lua_state, opts, 'complex')

            new_basis = new .or. aot_exists(lua_state, opts, 'int_file') &
                            .or. aot_exists(lua_state, opts, 'CAS')

            if (new_basis) then
                ! [todo] - deallocate existing basis info and start afresh.

                call init_system(sys)
                if (parent) call check_sys(sys)
                call read_in_integrals(sys, verbose=verbose)
                call init_generic_system_basis(sys)

                if (sys%momentum_space) then
                    call print_mom_sym_info(sys)
                else
                    call print_pg_symmetry_info(sys)
                end if
            end if
        end if

        if (parent) then
            ! This check can't go with the rest of check_sys as it has to happen after the basis is initialised.
            if (sys%symmetry < huge(0) .and. (sys%symmetry > sys%sym_max .or. sys%symmetry < sys%sym0)) &
                call stop_all('lua_read_in', 'Symmetry out of range.')
        end if

        call warn_unused_args(lua_state, keys, opts)

        call push_sys(lua_state, sys, new_basis)
        nreturn = 1

        call cpu_time(t2)
        call register_timing(lua_state, "Generic system initialisation", t2-t1)

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
        !        finite = true/false,
        !        triangular = true/false,
        !        staggered_magnetic_field = field,
        !        magnetic_field = field,
        !    }

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int
        use flu_binding, only: flu_State, flu_copyptr, flu_gettop
        use aot_top_module, only: aot_top_get_val
        use aot_table_module, only: aot_table_top, aot_get_val, aot_exists, aot_table_close

        use lua_hande_utils, only: warn_unused_args, register_timing
        use system, only: sys_t, heisenberg, init_system
        use basis, only: init_model_basis_fns
        use real_lattice, only: init_real_space
        use check_input, only: check_sys
        use parallel, only: parent

        integer(c_int) :: nreturn
        type(c_ptr), value :: L
        type(flu_State) :: lua_state

        type(sys_t), pointer :: sys
        integer :: opts
        real :: t1, t2
        logical :: new, new_basis, verbose
        integer :: err

        character(24), parameter :: keys(9) = [character(24) :: 'sys', 'ms', 'J', 'lattice', 'magnetic_field', &
                                                                'staggered_magnetic_field', 'triangular', 'finite', 'verbose']

        call cpu_time(t1)

        lua_state = flu_copyptr(L)
        call get_sys_t(lua_state, sys, new)

        ! get a handle to the table...
        opts = aot_table_top(lua_state)

        sys%system = heisenberg

        call set_common_sys_options(lua_state, sys, verbose, opts)
        call aot_get_val(sys%heisenberg%J, err, lua_state, opts, 'J')
        call aot_get_val(sys%heisenberg%magnetic_field, err, lua_state, opts, 'magnetic_field')
        call aot_get_val(sys%heisenberg%staggered_magnetic_field, err, lua_state, opts, 'staggered_magnetic_field')
        call aot_get_val(sys%real_lattice%finite_cluster, err, lua_state, opts, 'finite')
        call aot_get_val(sys%lattice%triangular_lattice, err, lua_state, opts, 'triangular')

        new_basis = aot_exists(lua_state, opts, 'lattice') .or. new

        if (new_basis) then
            call get_lattice(lua_state, sys, opts)
            ! [todo] - deallocate existing basis info and start afresh.
            call init_system(sys)
            if (parent) call check_sys(sys)
            call init_model_basis_fns(sys, verbose=verbose)
            call init_generic_system_basis(sys)
            call init_real_space(sys)
        end if

        call warn_unused_args(lua_state, keys, opts)

        call push_sys(lua_state, sys, new_basis)
        nreturn = 1

        call cpu_time(t2)
        call register_timing(lua_state, "Heisenberg model initialisation", t2-t1)

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
        !        twist = {...},    -- D-dimensional vector
        !    }
        !    Returns: sys_t object.

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int
        use flu_binding, only: flu_State, flu_copyptr, flu_gettop
        use aot_top_module, only: aot_top_get_val
        use aot_table_module, only: aot_table_top, aot_get_val, aot_exists, aot_table_close

        use lua_hande_utils, only: warn_unused_args, register_timing
        use system, only: sys_t, ueg, init_system
        use basis, only: init_model_basis_fns
        use momentum_symmetry, only: init_momentum_symmetry
        use ueg_system, only: init_ueg_proc_pointers
        use check_input, only: check_sys
        use parallel, only: parent
        use errors, only: stop_all

        integer(c_int) :: nreturn
        type(c_ptr), value :: L
        type(flu_State) :: lua_state

        type(sys_t), pointer :: sys
        integer :: opts
        real :: t1, t2
        logical :: new_basis, new, verbose
        integer :: err

        character(10), parameter :: keys(10) = [character(10) :: 'sys', 'cutoff', 'dim', 'rs', 'nel', 'electrons', &
                                               'ms', 'sym', 'twist', 'verbose']

        call cpu_time(t1)

        lua_state = flu_copyptr(L)
        call get_sys_t(lua_state, sys, new)

        ! Get a handle to the table...
        opts = aot_table_top(lua_state)

        sys%system = ueg
        sys%lattice%ndim = 3

        ! Parse table for options...
        call set_common_sys_options(lua_state, sys, verbose, opts)

        new_basis = aot_exists(lua_state, opts, 'cutoff') .or. &
                    aot_exists(lua_state, opts, 'dim')    .or. &
                    aot_exists(lua_state, opts, 'rs')    .or. &
                    aot_exists(lua_state, opts, 'nel')    .or. &
                    aot_exists(lua_state, opts, 'electrons') .or. new

        call aot_get_val(sys%ueg%ecutoff, err, lua_state, opts, 'cutoff')
        call aot_get_val(sys%ueg%r_s, err, lua_state, opts, 'rs')
        call aot_get_val(sys%lattice%ndim, err, lua_state, opts, 'dim')

        call get_ktwist(lua_state, sys, opts)

        if (new_basis) then
            ! [todo] - deallocate existing basis info and start afresh.

            call init_system(sys)
            if (parent) call check_sys(sys)
            call init_model_basis_fns(sys, verbose=verbose)
            call init_generic_system_basis(sys)
            call init_momentum_symmetry(sys)
            call init_ueg_proc_pointers(sys%lattice%ndim, sys%ueg)
        end if

        if (parent) then
            ! This check can't go with the rest of check_sys as it has to happen after the basis is initialised.
            if (sys%symmetry < huge(0) .and. (sys%symmetry > sys%sym_max .or. sys%symmetry < sys%sym0)) &
                call stop_all('lua_ueg', 'Symmetry out of range.')
        end if

        call warn_unused_args(lua_state, keys, opts)

        call push_sys(lua_state, sys, new_basis)
        nreturn = 1

        call cpu_time(t2)
        call register_timing(lua_state, "UEG initialisation", t2-t1)

    end function lua_ueg

    function lua_ringium(L) result(nreturn) bind(c)

        ! Create/modify a ringium system

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        !    ringium{
        !            sys = sys_old, -- new sys_t object is created if not passed
        !            electrons = N,
        !            radius = R
        !            maxlz = lzcutoff
        !    }
        !    Returns: sys_t object

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int
        use flu_binding, only: flu_state, flu_copyptr, flu_gettop
        use aot_top_module, only: aot_top_get_val
        use aot_table_module, only: aot_table_top, aot_get_val, aot_exists, aot_table_close

        use lua_hande_utils, only: warn_unused_args, register_timing
        use system, only: sys_t, ringium, init_system
        use basis, only: init_model_basis_fns
        use ringium_system, only: init_symmetry_ringium
        use check_input, only: check_sys
        use parallel, only: parent

        integer(c_int) :: nreturn
        type(c_ptr), value :: L
        type(flu_state) :: lua_state

        type(sys_t), pointer :: sys
        integer :: opts
        real :: t1, t2
        logical :: new_basis, new, verbose
        integer :: err

        character(10), parameter :: keys(7) = [character(10) :: 'sys', 'nel', 'electrons', 'radius', 'maxlz', 'sym', 'verbose']

        call cpu_time(t1)

        lua_state = flu_copyptr(L)
        call get_sys_t(lua_state, sys, new)

        ! Get a handle to the table...
        opts = aot_table_top(lua_state)

        sys%system = ringium
        sys%lattice%ndim = 1

        ! Parse table for options...
        call set_common_sys_options(lua_state, sys, verbose, opts)
        ! Enforce spin polarisation
        sys%ms = sys%nel

        new_basis = aot_exists(lua_state, opts, 'maxlz') .or. &
                    aot_exists(lua_state, opts, 'radius') .or. &
                    aot_exists(lua_state, opts, 'nel') .or. &
                    aot_exists(lua_state, opts, 'electrons') .or. new

        call aot_get_val(sys%ringium%radius, err, lua_state, opts, 'radius')
        call aot_get_val(sys%ringium%maxlz, err, lua_state, opts, 'maxlz')

        if (new_basis) then
            ! [todo] - deallocate existing basis info and start afresh.

            call init_system(sys)
            if (parent) call check_sys(sys)
            call init_model_basis_fns(sys, verbose=verbose)
            call init_generic_system_basis(sys)
            call init_symmetry_ringium(sys)
        end if

        call warn_unused_args(lua_state, keys, opts)

        call push_sys(lua_state, sys, new_basis)
        nreturn = 1

        call cpu_time(t2)
        call register_timing(lua_state, "Ringium initialisation", t2-t1)

    end function lua_ringium

    function lua_dealloc_sys(L) result(nresult) bind(c)

        ! Deallocate a sys object.  Expects to be called from lua with a single argument --
        ! the sys object to be deallocated.

        ! In/Out:
        !   L: lua state (bare C pointer).

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_f_pointer, c_null_ptr
        use flu_binding, only: flu_State, flu_copyptr, flu_pushstring, flu_pushlightuserdata, flu_settable
        use aot_table_ops_module, only: aot_table_top
        use aot_table_module, only: aot_get_val, aot_table_close

        use system, only: sys_t
        use dealloc, only: dealloc_sys_t
        use lua_hande_utils, only: get_userdata

        integer(c_int) :: nresult
        type(c_ptr), value :: L

        type(flu_State) :: lua_state
        integer :: sys_table
        type(c_ptr) :: sys_ptr
        type(sys_t), pointer :: sys

        lua_state = flu_copyptr(L)

        sys_table = aot_table_top(lua_state)
        call get_userdata(lua_state, sys_table, "sys", sys_ptr)
        call c_f_pointer(sys_ptr, sys)

        if (associated(sys)) then
            call dealloc_sys_t(sys)
            deallocate(sys)
        end if

        ! Update table with deallocated pointer.
        call flu_pushstring(lua_state, "sys")
        call flu_pushlightuserdata(lua_state, c_null_ptr)
        call flu_settable(lua_state, sys_table)

        call aot_table_close(lua_state, sys_table)

        nresult = 0

    end function lua_dealloc_sys

end module lua_hande_system
