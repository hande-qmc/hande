module lua_hande_calc

! Lua wrappers to for calculation procedures.  See top-level comments in 
! lua_hande for more details about working with the Lua API.

implicit none

contains

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
        type(c_ptr), value :: l

        type(flu_state) :: lua_state

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

    function lua_simple_fciqmc(L) result(nresult) bind(c)

        ! Run an FCIQMC calculation using the simple algorithm (slow, serial,
        ! memory hungry).

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        !    simple_fciqmc {
        !       sys = system,
        !       qmc = { ... },
        !       restart = { ... },
        !       reference = { ... },
        !    }

        ! See interface documentation for the relevant read_TYPE procedure to
        ! understand the options available within a given subtable.

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_f_pointer
        use flu_binding, only: flu_State, flu_copyptr
        use aot_table_module, only: aot_table_top, aot_table_close

        use simple_fciqmc, only: do_simple_fciqmc
        use lua_hande_system, only: get_sys_t
        use lua_hande_utils, only: warn_unused_args
        use qmc_data, only: qmc_in_t, restart_in_t, reference_t
        use qmc, only: init_proc_pointers
        use system, only: sys_t

        use calc, only: ms_in
        use system, only: set_spin_polarisation

        integer :: nresult
        type(c_ptr), value :: L

        type(flu_state) :: lua_state
        type(sys_t), pointer :: sys
        type(qmc_in_t) :: qmc_in
        type(restart_in_t) :: restart_in
        type(reference_t) :: reference

        integer :: opts
        character(10), parameter :: keys(4) = [character(10) :: 'sys', 'qmc', 'restart', 'reference']

        lua_state = flu_copyptr(L)
        call get_sys_t(lua_state, sys)
        ! [todo] - do spin polarisation in system setup.
        call set_spin_polarisation(sys%basis%nbasis, ms_in, sys)

        ! Get main table.
        opts = aot_table_top(lua_state)
        call read_qmc_in(lua_state, opts, qmc_in)
        ! [todo] - implement
        !call read_restart_in(lua_state, opts, restart_in)
        call read_reference_t(lua_state, opts, sys, reference)
        reference%ex_level = sys%nel
        call warn_unused_args(lua_state, keys, opts)
        call aot_table_close(lua_state, opts)

        call do_simple_fciqmc(sys, qmc_in, restart_in, reference)

        nresult = 0

    end function lua_simple_fciqmc

    function lua_fciqmc(L) result(nresult) bind(c)

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        !    fciqmc {
        !       sys = system,
        !       qmc = { ... },
        !       fciqmc = { ... },
        !       semi_stoch = { ... },
        !       restart = { ... },
        !       load_bal = { ... },
        !       reference = { ... },
        !    }

        ! See interface documentation for the relevant read_TYPE procedure to
        ! understand the options available within a given subtable.

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_f_pointer
        use flu_binding, only: flu_State, flu_copyptr
        use aot_table_module, only: aot_table_top, aot_table_close

        use dmqmc_data, only: dmqmc_in_t
        use fciqmc, only: do_fciqmc
        use lua_hande_system, only: get_sys_t
        use lua_hande_utils, only: warn_unused_args
        use qmc_data, only: qmc_in_t, fciqmc_in_t, semi_stoch_in_t, restart_in_t, load_bal_in_t, reference_t
        use qmc, only: init_proc_pointers
        use system, only: sys_t

        use calc, only: ms_in
        use system, only: set_spin_polarisation

        integer :: nresult
        type(c_ptr), value :: L

        type(flu_state) :: lua_state
        type(sys_t), pointer :: sys
        type(qmc_in_t) :: qmc_in
        type(fciqmc_in_t) :: fciqmc_in
        type(semi_stoch_in_t) :: semi_stoch_in
        type(restart_in_t) :: restart_in
        type(load_bal_in_t) :: load_bal_in
        type(reference_t) :: reference

        type(dmqmc_in_t) :: dmqmc_defaults
        integer :: opts
        character(10), parameter :: keys(7) = [character(10) :: 'sys', 'qmc', 'fciqmc', 'semi_stoch', 'restart', &
                                                                'load_bal', 'reference']

        lua_state = flu_copyptr(L)
        call get_sys_t(lua_state, sys)
        ! [todo] - do spin polarisation in system setup.
        call set_spin_polarisation(sys%basis%nbasis, ms_in, sys)

        ! Get main table.
        opts = aot_table_top(lua_state)

        call read_qmc_in(lua_state, opts, qmc_in)
        call read_fciqmc_in(lua_state, opts, fciqmc_in)
        call read_semi_stoch_in(lua_state, opts, semi_stoch_in)
        ! [todo] - implement
        !call read_restart_in(lua_state, opts, restart_in)
        call read_load_bal_in(lua_state, opts, load_bal_in)
        call read_reference_t(lua_state, opts, sys, reference)
        reference%ex_level = sys%nel
        call warn_unused_args(lua_state, keys, opts)
        call aot_table_close(lua_state, opts)

        call init_proc_pointers(sys, qmc_in, dmqmc_defaults, reference)
        call do_fciqmc(sys, qmc_in, fciqmc_in, semi_stoch_in, restart_in, load_bal_in, reference)

        nresult = 0

    end function lua_fciqmc

    function lua_ccmc(L) result(nresult) bind(c)

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        !    ccmc{
        !       sys = system,
        !       qmc = { ... },
        !       ccmc = { ... }
        !       restart = { ... },
        !       reference = { ... },
        !    }

        ! See interface documentation for the relevant read_TYPE procedure to
        ! understand the options available within a given subtable.

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_f_pointer
        use flu_binding, only: flu_State, flu_copyptr
        use aot_table_module, only: aot_table_top, aot_table_close

        use dmqmc_data, only: dmqmc_in_t
        use ccmc, only: do_ccmc
        use lua_hande_system, only: get_sys_t
        use lua_hande_utils, only: warn_unused_args
        use qmc_data, only: qmc_in_t, ccmc_in_t, semi_stoch_in_t, restart_in_t, load_bal_in_t, reference_t
        use qmc, only: init_proc_pointers
        use system, only: sys_t

        use calc, only: ms_in
        use system, only: set_spin_polarisation

        integer :: nresult
        type(c_ptr), value :: L

        type(flu_state) :: lua_state
        type(sys_t), pointer :: sys
        type(qmc_in_t) :: qmc_in
        type(ccmc_in_t) :: ccmc_in
        type(semi_stoch_in_t) :: semi_stoch_in
        type(restart_in_t) :: restart_in
        type(load_bal_in_t) :: load_bal_in
        type(reference_t) :: reference

        type(dmqmc_in_t) :: dmqmc_defaults
        integer :: opts
        character(10), parameter :: keys(5) = [character(10) :: 'sys', 'qmc', 'ccmc', 'restart', 'reference']

        lua_state = flu_copyptr(L)
        call get_sys_t(lua_state, sys)
        ! [todo] - do spin polarisation in system setup.
        call set_spin_polarisation(sys%basis%nbasis, ms_in, sys)

        ! Get main table.
        opts = aot_table_top(lua_state)
        call read_qmc_in(lua_state, opts, qmc_in)
        call read_ccmc_in(lua_state, opts, ccmc_in)
        ! note that semi-stochastic is not (yet) available in CCMC.
        !call read_semi_stoch_in(lua_state, opts, semi_stoch_in)
        ! [todo] - implement
        !call read_restart_in(lua_state, opts, restart_in)
        ! load balancing is not available in CCMC; must use default settings.
        call read_reference_t(lua_state, opts, sys, reference)
        reference%ex_level = sys%nel
        call warn_unused_args(lua_state, keys, opts)
        call aot_table_close(lua_state, opts)

        call init_proc_pointers(sys, qmc_in, dmqmc_defaults, reference)
        call do_ccmc(sys, qmc_in, ccmc_in, semi_stoch_in, restart_in, load_bal_in, reference)

        nresult = 0

    end function lua_ccmc

    ! --- table->derived type wrappers ---

    subroutine read_qmc_in(lua_state, opts, qmc_in)

        ! Read in a qmc table to a qmc_in_t object.

        ! qmc = {
        !     tau = tau,                                  -- required
        !     init_pop = N,                               -- required
        !     mc_cycles = ncycles,                        -- required
        !     nreports = nreport,                         -- required
        !     state_size = walker_length,                 -- required
        !     spawned_state_size = spawned_walker_length, -- required
        !     target_population = pop,
        !     real_amplitudes = true/false,
        !     spawn_cutoff = cutoff,
        !     no_renorm = true/false,
        !     tau_search = true/false,
        !     pattempt_single = prob,
        !     pattempt_double = prob,
        !     initial_shift = shift,
        !     shift_damping = damp_factor,
        !     initiator = true/false,
        !     initiator_threshold = pop,
        ! }

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !    opts: handle for the table containing the qmc table.
        ! Out:
        !    qmc_in: qmc_in_t object containing generic QMC input options.

        use flu_binding, only: flu_State
        use aot_table_module, only: aot_get_val, aot_exists, aot_table_open, aot_table_close

        use calc, only: GLOBAL_META, gen_seed
        use qmc_data, only: qmc_in_t
        use lua_hande_utils, only: warn_unused_args
        use errors, only: stop_all

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        type(qmc_in_t), intent(out) :: qmc_in

        integer :: qmc_table, err
        character(len=10) :: str

        call aot_table_open(lua_state, opts, qmc_table, 'qmc')

        ! Required arguments
        call aot_get_val(qmc_in%tau, err, lua_state, qmc_table, 'tau')
        if (err /= 0) call stop_all('read_qmc_in', 'tau not set.')
        call aot_get_val(qmc_in%D0_population, err, lua_state, qmc_table, 'init_pop')
        if (err /= 0) call stop_all('read_qmc_in', 'init_pop not set.')
        call aot_get_val(qmc_in%ncycles, err, lua_state, qmc_table, 'mc_cycles')
        if (err /= 0) call stop_all('read_qmc_in', 'mc_cycles not set.')
        call aot_get_val(qmc_in%nreport, err, lua_state, qmc_table, 'nreports')
        if (err /= 0) call stop_all('read_qmc_in', 'nreports not set.')
        call aot_get_val(qmc_in%walker_length, err, lua_state, qmc_table, 'state_size')
        if (err /= 0) call stop_all('read_qmc_in', 'state_size not set.')
        call aot_get_val(qmc_in%spawned_walker_length, err, lua_state, qmc_table, 'spawned_state_size')
        if (err /= 0) call stop_all('read_qmc_in', 'spawned_state_size not set.')

        ! Optional arguments (defaults set in derived type).
        call aot_get_val(qmc_in%real_amplitudes, err, lua_state, qmc_table, 'real_amplitudes')
        call aot_get_val(qmc_in%spawn_cutoff, err, lua_state, qmc_table, 'spawn_cutoff')
        call aot_get_val(qmc_in%no_renorm, err, lua_state, qmc_table, 'no_renorm')
        call aot_get_val(qmc_in%pattempt_single, err, lua_state, qmc_table, 'pattempt_single')
        call aot_get_val(qmc_in%pattempt_double, err, lua_state, qmc_table, 'pattempt_double')
        call aot_get_val(qmc_in%tau_search, err, lua_state, qmc_table, 'tau_search')
        call aot_get_val(qmc_in%initial_shift, err, lua_state, qmc_table, 'initial_shift')
        call aot_get_val(qmc_in%shift_damping, err, lua_state, qmc_table, 'shift_damping')
        call aot_get_val(qmc_in%target_particles, err, lua_state, qmc_table, 'target_population')
        call aot_get_val(qmc_in%initiator_approx, err, lua_state, qmc_table, 'initiator')
        call aot_get_val(qmc_in%initiator_pop, err, lua_state, qmc_table, 'initiator_threshold')

        ! Optional arguments requiring special care.
        if (aot_exists(lua_state, qmc_table, 'rng_seed')) then
            call aot_get_val(qmc_in%seed, err, lua_state, qmc_table, 'rng_seed')
        else
            qmc_in%seed = gen_seed(GLOBAL_META%uuid)
        end if
        if (aot_exists(lua_state, qmc_table, 'vary_shift_from')) then
            call aot_get_val(qmc_in%vary_shift_from, err, lua_state, qmc_table, 'vary_shift_from')
            if (err /= 0) then
                ! Perhaps passed in a string?
                call aot_get_val(str, err, lua_state, qmc_table, 'vary_shift_from')
                qmc_in% vary_shift_from_proje = trim(str) == 'proje'
            end if
        end if

        ! [todo] - check unused args.
        call warn_unused_args(lua_state, ['tau'], qmc_table)

        call aot_table_close(lua_state, qmc_table)

    end subroutine read_qmc_in

    subroutine read_fciqmc_in(lua_state, opts, fciqmc_in)

        ! Read in an fciqmc table (if it exists) to an fciqmc_in object.

        ! fciqmc = {
        !     non_blocking_comm = true/false,
        !     load_balancing = true/false,
        !     init_spin_inverse_reference_det = true/false,
        !     select_reference_det = true, -- OR
        !     select_reference_det = {
        !        update_every = mc_cycles,
        !        pop_factor = factor,
        !     }
        ! }

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !    opts: handle for the table containing the fciqmc table.
        ! Out:
        !    fciqmc_in: fciqmc_in_t object containing FCIQMC-specific input options.

        use flu_binding, only: flu_State
        use aot_table_module, only: aot_get_val, aot_exists, aot_table_open, aot_table_close

        use qmc_data, only: fciqmc_in_t
        use lua_hande_utils, only: warn_unused_args

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        type(fciqmc_in_t), intent(out) :: fciqmc_in

        integer :: fciqmc_table, ref_det, err
        character(len=10) :: str
        logical :: ref_det_flag

        if (aot_exists(lua_state, opts, 'fciqmc')) then

            call aot_table_open(lua_state, opts, fciqmc_table, 'fciqmc')

            ! Optional arguments (defaults set in derived type).
            call aot_get_val(fciqmc_in%non_blocking_comm, err, lua_state, fciqmc_table, 'non_blocking_comm')
            call aot_get_val(fciqmc_in%doing_load_balancing, err, lua_state, fciqmc_table, 'load_balancing')
            call aot_get_val(fciqmc_in%init_spin_inv_D0, err, lua_state, fciqmc_table, 'init_spin_inverse_reference_det')

            ! Optional arguments requiring special care.
            if (aot_exists(lua_state, fciqmc_table, 'select_reference_det')) then
                call aot_table_open(lua_state, fciqmc_table, ref_det, 'select_reference_det')
                if (ref_det == 0) then
                    ! Just passed a boolean (hopefully).
                    call aot_get_val(ref_det_flag, err, lua_state, ref_det, 'select_reference_det', default=.false.)
                    if (ref_det_flag) fciqmc_in%select_ref_det_every_nreports = 20
                else
                    call aot_table_close(lua_state, ref_det)
                    call aot_get_val(fciqmc_in%select_ref_det_every_nreports, err, lua_state, ref_det, 'update_every', default=20)
                    call aot_get_val(fciqmc_in%ref_det_factor, err, lua_state, ref_det, 'pop_factor')
                end if
            end if

            ! [todo] - check unused args.
            !call warn_unused_args(lua_state, [], fciqmc_table)

            call aot_table_close(lua_state, fciqmc_table)

        end if

    end subroutine read_fciqmc_in

    subroutine read_semi_stoch_in(lua_state, opts, semi_stoch_in)

        ! Read in an semi_stoch table (if it exists) to an semi_stoch_in object.

        ! semi_stoch = {
        !     space = 'read'/'high',            -- required if semi_stoch table exists.
        !     size = target_pop,                -- required if space == 'high'
        !     start_iteration = mc_cycles,
        !     shift_start_iteration = mc_cycle, -- overrides start_iteration
        !     separate_annihilation = true/false,
        !     write_determ_space = true/false,
        ! }

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !    opts: handle for the table containing the semi_stoch table.
        ! Out:
        !    semi_stoch_in: semi_stoch_in_t object containing semi-stochastic-specific input options.

        use flu_binding, only: flu_State
        use aot_table_module, only: aot_get_val, aot_exists, aot_table_open, aot_table_close
        use errors, only: stop_all

        use qmc_data, only: semi_stoch_in_t, high_pop_determ_space, read_determ_space
        use lua_hande_utils, only: warn_unused_args

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        type(semi_stoch_in_t), intent(out) :: semi_stoch_in

        integer :: semi_stoch_table, ref_det, err
        character(len=10) :: str
        logical :: ref_det_flag

        if (aot_exists(lua_state, opts, 'semi_stoch')) then

            call aot_table_open(lua_state, opts, semi_stoch_table, 'semi_stoch')

            ! Required arguments.
            if (aot_exists(lua_state, semi_stoch_table, 'space')) then
                call aot_get_val(str, err, lua_state, semi_stoch_table, 'space')
                select case(str)
                case('read')
                    semi_stoch_in%space_type = read_determ_space
                    semi_stoch_in%target_size = -1
                case('high')
                    semi_stoch_in%space_type = high_pop_determ_space
                    call aot_get_val(semi_stoch_in%target_size, err, lua_state, semi_stoch_table, 'size')
                    if (err /= 0) call stop_all('read_semi_stoch_in', 'Target space size not given.')
                end select
            else
                call stop_all('read_semi_stoch_in', 'Deterministic space not specified.')
            end if

            ! Optional arguments (defaults set in derived type).
            call aot_get_val(semi_stoch_in%start_iter, err, lua_state, semi_stoch_table, 'start_iteration')
            call aot_get_val(semi_stoch_in%write_determ_space, err, lua_state, semi_stoch_table, 'write_determ_space')
            call aot_get_val(semi_stoch_in%separate_annihil, err, lua_state, semi_stoch_table, 'separate_annihilation')

            ! Optional arguments requiring special care.
            if (aot_exists(lua_state, semi_stoch_table, 'shift_start_iteration')) then
                semi_stoch_in%start_iter = huge(0) ! fixed up once the shift comes on...
                call aot_get_val(semi_stoch_in%shift_iter, err, lua_state, semi_stoch_table, 'shift_start_iteration')
            end if

            ! [todo] - check unused args.
            !call warn_unused_args(lua_state, [], semi_stoch_table)

            call aot_table_close(lua_state, semi_stoch_table)

        end if

    end subroutine read_semi_stoch_in

    subroutine read_ccmc_in(lua_state, opts, ccmc_in)

        ! Read in an ccmc table (if it exists) to an ccmc_in object.

        ! ccmc = {
        !     move_freq = frequency,
        !     cluster_multispawn_threhold = threshold,
        !     full_non_composite = true/false,
        !     linked = true/false,
        ! }

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !    opts: handle for the table containing the ccmc table.
        ! Out:
        !    ccmc_in: ccmc_in_t object containing ccmc-specific input options.

        use flu_binding, only: flu_State
        use aot_table_module, only: aot_get_val, aot_exists, aot_table_open, aot_table_close

        use qmc_data, only: ccmc_in_t
        use lua_hande_utils, only: warn_unused_args

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        type(ccmc_in_t), intent(out) :: ccmc_in

        integer :: ccmc_table, ref_det, err
        character(len=10) :: str
        logical :: ref_det_flag

        if (aot_exists(lua_state, opts, 'ccmc')) then

            call aot_table_open(lua_state, opts, ccmc_table, 'ccmc')

            ! Optional arguments (defaults set in derived type).
            call aot_get_val(ccmc_in%move_freq, err, lua_state, ccmc_table, 'move_frequency')
            call aot_get_val(ccmc_in%cluster_multispawn_threshold, err, lua_state, ccmc_table, 'cluster_multispawn_threhold')
            call aot_get_val(ccmc_in%full_nc, err, lua_state, ccmc_table, 'full_non_composite')
            call aot_get_val(ccmc_in%linked, err, lua_state, ccmc_table, 'linked')

            ! [todo] - check unused args.
            !call warn_unused_args(lua_state, [], ccmc_table)

            call aot_table_close(lua_state, ccmc_table)

        end if

    end subroutine read_ccmc_in

    subroutine read_load_bal_in(lua_state, opts, load_bal_in)

        ! Read in an load_bal table (if it exists) to an load_bal_in object.

        ! load_bal = {
        !     nslots = N,
        !     min_pop = pop,
        !     target = percent,
        !     max_attempts = max,
        !     write = true/false,
        ! }

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !    opts: handle for the table containing the load_bal table.
        ! Out:
        !    load_bal_in: load_bal_in_t object containing load-balancing-specific input options.

        use flu_binding, only: flu_State
        use aot_table_module, only: aot_get_val, aot_exists, aot_table_open, aot_table_close

        use qmc_data, only: load_bal_in_t
        use lua_hande_utils, only: warn_unused_args

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        type(load_bal_in_t), intent(out) :: load_bal_in

        integer :: load_bal_table, ref_det, err
        character(len=10) :: str
        logical :: ref_det_flag

        if (aot_exists(lua_state, opts, 'load_bal')) then

            call aot_table_open(lua_state, opts, load_bal_table, 'load_bal')

            ! Optional arguments (defaults set in derived type).
            call aot_get_val(load_bal_in%nslots, err, lua_state, load_bal_table, 'nslots')
            call aot_get_val(load_bal_in%pop, err, lua_state, load_bal_table, 'min_pop')
            call aot_get_val(load_bal_in%percent, err, lua_state, load_bal_table, 'target')
            call aot_get_val(load_bal_in%max_attempts, err, lua_state, load_bal_table, 'max_attempts')
            call aot_get_val(load_bal_in%write_info, err, lua_state, load_bal_table, 'write')

            ! [todo] - check unused args.
            !call warn_unused_args(lua_state, [], load_bal_table)

            call aot_table_close(lua_state, load_bal_table)

        end if

    end subroutine read_load_bal_in

    subroutine read_reference_t(lua_state, opts, sys, ref)

        ! Read in a reference table (if it exists) to a reference_t object.

        ! reference = {
        !     det = { ... }, -- N-electron vector
        !     hilbert_space_det = { ... }, -- N-electron vector
        !     ex_level = truncation_level,
        ! }

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !    opts: handle for the table containing the reference table.
        !    sys: sys_t object describing the system of interest.
        ! Out:
        !    ref: reference_t object contianing the input options describing the
        !        reference.  Note that ex_level is set to the number of electrons
        !        if not provided (incl. if the reference table is not present in
        !        opts).

        use flu_binding, only: flu_State
        use aot_vector_module, only: aot_get_val
        use aot_table_module, only: aot_get_val, aot_exists, aot_table_open, aot_table_close

        use qmc_data, only: reference_t
        use system, only: sys_t

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        type(sys_t), intent(in) :: sys
        type(reference_t), intent(out) :: ref

        integer :: ref_table, err
        integer, allocatable :: err_arr(:)

        ! Set to full space by default.
        ref%ex_level = sys%nel

        if (aot_exists(lua_state, opts, 'reference')) then

            call aot_table_open(lua_state, opts, ref_table, 'reference')

            ! Optional arguments requiring special handling.

            ! Check if deeterminants exist first to avoid aotus allocating components unnecessarily.
            if (aot_exists(lua_state, ref_table, 'det')) then
                call aot_get_val(ref%occ_list0, err_arr, sys%nel, lua_state, ref_table, 'det')
            end if
            if (aot_exists(lua_state, ref_table, 'hilbert_space_det')) then
                call aot_get_val(ref%hs_occ_list0, err_arr, sys%nel, lua_state, ref_table, 'hilbert_space_det')
            end if

            call aot_get_val(ref%ex_level, err, lua_state, ref_table, 'ex_level')

            ! [todo] - check unused args.
            !call warn_unused_args(lua_state, [], load_bal_table)

            call aot_table_close(lua_state, ref_table)

        end if

    end subroutine read_reference_t

end module lua_hande_calc
