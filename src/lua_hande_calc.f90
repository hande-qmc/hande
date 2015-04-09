module lua_hande_calc

! Lua wrappers to for calculation procedures.  See top-level comments in 
! lua_hande for more details about working with the Lua API.

implicit none

contains

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
        use lua_hande_system, only: get_sys_t
        use lua_hande_utils, only: warn_unused_args
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
        use lua_hande_system, only: get_sys_t
        use lua_hande_utils, only: warn_unused_args
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

end module lua_hande_calc
