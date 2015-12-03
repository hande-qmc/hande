module lua_hande_calc

! Lua wrappers to for calculation procedures.  See top-level comments in
! lua_hande for more details about working with the Lua API.

implicit none

contains

    ! --- Calculation wrappers ---

    function lua_fci(L) result(nresult) bind(c)

        ! Run an FCI calculation.

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        !    fci {
        !           sys = system,        -- required
        !           fci = { ... },
        !           lanczos = { ... },
        !           reference = { ... },
        !    }

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int
        use flu_binding, only: flu_State, flu_copyptr
        use aot_table_module, only: aot_table_top, aot_exists, aot_table_close

        use calc, only: calc_type, exact_diag, lanczos_diag
        use fci_lanczos, only: do_fci_lanczos
        use fci_lapack, only: do_fci_lapack
        use fci_utils, only: fci_in_t
        use lua_hande_system, only: get_sys_t
        use lua_hande_utils, only: warn_unused_args
        use qmc_data, only: reference_t
        use system, only: sys_t

        integer(c_int) :: nresult
        type(c_ptr), value :: L

        type(flu_State) :: lua_state
        integer :: opts
        type(sys_t), pointer :: sys
        type(fci_in_t) :: fci_in
        type(reference_t) :: ref
        logical :: use_sparse_hamil, lanczos
        character(12), parameter :: keys(4) = [character(12) :: 'sys', 'fci', 'lanczos', 'reference']

        lua_state = flu_copyptr(l)
        call get_sys_t(lua_state, sys)

        opts = aot_table_top(lua_state)
        lanczos = aot_exists(lua_state, opts, 'lanczos')
        call read_fci_in(lua_state, opts, sys%basis, fci_in, use_sparse_hamil)
        call read_reference_t(lua_state, opts, sys, ref)
        call warn_unused_args(lua_state, keys, opts)
        call aot_table_close(lua_state, opts)

        if (lanczos) then
            calc_type = lanczos_diag
            call do_fci_lanczos(sys, fci_in, ref, use_sparse_hamil)
        else
            calc_type = exact_diag
            call do_fci_lapack(sys, fci_in, ref)
        end if

        nresult = 0

    end function lua_fci

    function lua_hilbert_space(L) result(nresult) bind(c)

        ! Run a Monte Carlo calculation to estimate the size of the Hilbert
        ! space.

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        !    hilbert_space {
        !       sys = system,      -- required
        !       hilbert = { ... }, -- required
        !    }

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_f_pointer

        use flu_binding, only: flu_State, flu_copyptr
        use aot_table_module, only: aot_table_top, aot_exists, aot_table_close

        use parallel, only: parent
        use errors, only: stop_all
        use lua_hande_system, only: get_sys_t
        use lua_hande_utils, only: warn_unused_args
        use system, only: sys_t

        use calc, only: calc_type, mc_hilbert_space
        use hilbert_space, only: estimate_hilbert_space

        integer(c_int) :: nresult
        type(c_ptr), value :: L

        type(flu_State) :: lua_state

        type(sys_t), pointer :: sys
        integer :: truncation_level, nattempts, ncycles, rng_seed
        integer, allocatable :: ref_det(:)
        integer :: opts
        character(12), parameter :: keys(2) = [character(12) :: 'sys', 'hilbert']

        lua_state = flu_copyptr(l)
        call get_sys_t(lua_state, sys)

        opts = aot_table_top(lua_state)
        call read_hilbert_args(lua_state, opts, sys%nel, nattempts, ncycles, truncation_level, ref_det, rng_seed)
        call warn_unused_args(lua_state, keys, opts)
        call aot_table_close(lua_state, opts)

        ! AOTUS returns a vector of size 0 to denote a non-existent vector.
        if (size(ref_det) == 0) deallocate(ref_det)
        if (allocated(ref_det)) then
            if (size(ref_det) /= sys%nel .and. parent) call stop_all('lua_hilbert_space', &
                            'Reference determinant does not match the number of electrons in system.')
        end if

        calc_type = mc_hilbert_space
        call estimate_hilbert_space(sys, truncation_level, nattempts, ncycles, ref_det, rng_seed)

        ! [todo] - return estimate of space and error to lua.
        nresult = 0

    end function lua_hilbert_space

    function lua_canonical_energy(L) result(nresult) bind(c)

        ! Run a Monte Carlo calculation to estimate various mean-field energies of a
        ! system in the canonical ensemble.

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        !    canonical_energy {
        !       sys = system,      -- required
        !       canonical_energy = { ... }, -- required
        !    }

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_f_pointer
        use const, only: p

        use flu_binding, only: flu_State, flu_copyptr
        use aot_table_module, only: aot_table_top, aot_table_close

        use parallel, only: parent
        use errors, only: stop_all
        use lua_hande_system, only: get_sys_t
        use lua_hande_utils, only: warn_unused_args
        use system, only: sys_t

        use calc, only: calc_type, mc_canonical_energy_estimates
        use canonical_energy_estimates, only: estimate_canonical_energy

        integer(c_int) :: nresult
        type(c_ptr), value :: l

        type(flu_state) :: lua_state

        type(sys_t), pointer :: sys
        integer :: opts, rng_seed, ncycles, nattempts
        logical :: fermi_temperature, all_spin_sectors
        real(p) :: beta
        character(16), parameter :: keys(2) = [character(16) :: 'sys', 'canonical_energy']

        lua_state = flu_copyptr(L)
        call get_sys_t(lua_state, sys)

        if (.not. sys%chem_pot < huge(1.0_p) .and. parent) call stop_all('lua_canonical_energy', &
                                                                        'chem_pot: chemical potential not supplied.')

        opts = aot_table_top(lua_state)
        call read_canonical_energy_args(lua_state, opts, fermi_temperature, beta, nattempts, ncycles, rng_seed, all_spin_sectors)
        call warn_unused_args(lua_state, keys, opts)
        call aot_table_close(lua_state, opts)

        calc_type = mc_canonical_energy_estimates
        call estimate_canonical_energy(sys, fermi_temperature, beta, nattempts, ncycles, all_spin_sectors, rng_seed)

        ! [todo] - return estimate of various canonical mean-field energies and error to lua.
        nresult = 0

    end function lua_canonical_energy

    function lua_simple_fciqmc(L) result(nresult) bind(c)

        ! Run an FCIQMC calculation using the simple algorithm (slow, serial,
        ! memory hungry).

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        !    simple_fciqmc {
        !       sys = system,         -- required
        !       qmc = { ... },        -- required
        !       restart = { ... },
        !       reference = { ... },
        !       sparse = true/false,
        !    }

        ! See interface documentation for the relevant read_TYPE procedure to
        ! understand the options available within a given subtable.

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_f_pointer
        use flu_binding, only: flu_State, flu_copyptr
        use aot_table_module, only: aot_table_top, aot_table_close, aot_get_val

        use simple_fciqmc, only: do_simple_fciqmc
        use lua_hande_system, only: get_sys_t
        use lua_hande_utils, only: warn_unused_args
        use qmc_data, only: qmc_in_t, restart_in_t, reference_t
        use qmc, only: init_proc_pointers
        use system, only: sys_t

        use calc, only: calc_type, simple_fciqmc_calc, fciqmc_calc
        use system, only: set_spin_polarisation

        integer(c_int) :: nresult
        type(c_ptr), value :: L

        type(flu_state) :: lua_state
        type(sys_t), pointer :: sys
        type(qmc_in_t) :: qmc_in
        type(restart_in_t) :: restart_in
        type(reference_t) :: reference
        logical :: use_sparse_hamil

        integer :: opts, err
        character(12), parameter :: keys(5) = [character(12) :: 'sys', 'qmc', 'restart', 'reference', 'sparse']

        lua_state = flu_copyptr(L)
        call get_sys_t(lua_state, sys)
        ! [todo] - do spin polarisation in system setup.
        call set_spin_polarisation(sys%basis%nbasis, sys)

        ! Get main table.
        opts = aot_table_top(lua_state)
        call read_qmc_in(lua_state, opts, qmc_in, .true.)
        call read_restart_in(lua_state, opts, restart_in)
        call read_reference_t(lua_state, opts, sys, reference)
        call aot_get_val(use_sparse_hamil, err, lua_state, opts, 'sparse', default=.true.)
        call warn_unused_args(lua_state, keys, opts)
        call aot_table_close(lua_state, opts)

        calc_type = simple_fciqmc_calc + fciqmc_calc
        call do_simple_fciqmc(sys, qmc_in, restart_in, reference, use_sparse_hamil)

        nresult = 0

    end function lua_simple_fciqmc

    function lua_fciqmc(L) result(nresult) bind(c)

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        !    fciqmc {
        !       sys = system,          -- required
        !       qmc = { ... },         -- required
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

        use calc, only: calc_type, fciqmc_calc
        use system, only: set_spin_polarisation

        integer(c_int)  :: nresult
        type(c_ptr), value :: L

        type(flu_state) :: lua_state
        type(sys_t), pointer :: sys
        type(qmc_in_t) :: qmc_in
        type(fciqmc_in_t) :: fciqmc_in
        type(semi_stoch_in_t) :: semi_stoch_in
        type(restart_in_t) :: restart_in
        type(load_bal_in_t) :: load_bal_in
        type(reference_t) :: reference

        integer :: opts
        character(10), parameter :: keys(7) = [character(10) :: 'sys', 'qmc', 'fciqmc', 'semi_stoch', 'restart', &
                                                                'load_bal', 'reference']

        lua_state = flu_copyptr(L)
        call get_sys_t(lua_state, sys)
        ! [todo] - do spin polarisation in system setup.
        call set_spin_polarisation(sys%basis%nbasis, sys)

        ! Get main table.
        opts = aot_table_top(lua_state)

        call read_qmc_in(lua_state, opts, qmc_in)
        call read_fciqmc_in(lua_state, opts, fciqmc_in)
        call read_semi_stoch_in(lua_state, opts, qmc_in, semi_stoch_in)
        call read_restart_in(lua_state, opts, restart_in)
        call read_load_bal_in(lua_state, opts, load_bal_in)
        call read_reference_t(lua_state, opts, sys, reference)
        call warn_unused_args(lua_state, keys, opts)
        call aot_table_close(lua_state, opts)

        calc_type = fciqmc_calc
        call init_proc_pointers(sys, qmc_in, reference, fciqmc_in=fciqmc_in)
        call do_fciqmc(sys, qmc_in, fciqmc_in, semi_stoch_in, restart_in, load_bal_in, reference)

        nresult = 0

    end function lua_fciqmc

    function lua_ccmc(L) result(nresult) bind(c)

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        !    ccmc{
        !       sys = system,         --  required
        !       qmc = { ... },        --  required
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

        use calc, only: calc_type, ccmc_calc
        use system, only: set_spin_polarisation

        integer(c_int) :: nresult
        type(c_ptr), value :: L

        type(flu_state) :: lua_state
        type(sys_t), pointer :: sys
        type(qmc_in_t) :: qmc_in
        type(ccmc_in_t) :: ccmc_in
        type(semi_stoch_in_t) :: semi_stoch_in
        type(restart_in_t) :: restart_in
        type(load_bal_in_t) :: load_bal_in
        type(reference_t) :: reference

        integer :: opts
        character(10), parameter :: keys(5) = [character(10) :: 'sys', 'qmc', 'ccmc', 'restart', 'reference']

        lua_state = flu_copyptr(L)
        call get_sys_t(lua_state, sys)
        ! [todo] - do spin polarisation in system setup.
        call set_spin_polarisation(sys%basis%nbasis, sys)

        ! Get main table.
        opts = aot_table_top(lua_state)
        call read_qmc_in(lua_state, opts, qmc_in)
        call read_ccmc_in(lua_state, opts, ccmc_in)
        ! note that semi-stochastic is not (yet) available in CCMC.
        !call read_semi_stoch_in(lua_state, opts, qmc_in, semi_stoch_in)
        call read_restart_in(lua_state, opts, restart_in)
        ! load balancing is not available in CCMC; must use default settings.
        call read_reference_t(lua_state, opts, sys, reference)
        call warn_unused_args(lua_state, keys, opts)
        call aot_table_close(lua_state, opts)

        calc_type = ccmc_calc
        call init_proc_pointers(sys, qmc_in, reference)
        call do_ccmc(sys, qmc_in, ccmc_in, semi_stoch_in, restart_in, load_bal_in, reference)

        nresult = 0

    end function lua_ccmc

    function lua_dmqmc(L) result(nresult) bind(c)

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        !    dmqmc{
        !       sys = system,         -- required
        !       qmc = { ... },        -- required
        !       dmqmc = { ... },
        !       ipdmqmc = { ... },
        !       operators = { ... },
        !       rdm = { ... },
        !       restart = { ... },
        !       reference = { ... },
        !    }

        ! See interface documentation for the relevant read_TYPE procedure to
        ! understand the options available within a given subtable.

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_f_pointer
        use flu_binding, only: flu_State, flu_copyptr
        use aot_table_module, only: aot_table_top, aot_table_close

        use dmqmc_data, only: dmqmc_in_t, dmqmc_estimates_t
        use dmqmc, only: do_dmqmc
        use lua_hande_system, only: get_sys_t
        use lua_hande_utils, only: warn_unused_args
        use qmc_data, only: qmc_in_t, restart_in_t, load_bal_in_t, reference_t
        use qmc, only: init_proc_pointers
        use system, only: sys_t, heisenberg, read_in, ueg

        use calc, only: calc_type, dmqmc_calc, doing_dmqmc_calc, dmqmc_energy_squared
        use checking, only: check_allocate
        use real_lattice, only: create_next_nearest_orbs
        use system, only: set_spin_polarisation

        integer(c_int) :: nresult
        type(c_ptr), value :: L

        type(flu_state) :: lua_state
        type(sys_t), pointer :: sys
        type(qmc_in_t) :: qmc_in
        type(dmqmc_in_t) :: dmqmc_in
        type(dmqmc_estimates_t) :: dmqmc_estimates
        type(restart_in_t) :: restart_in
        type(load_bal_in_t) :: load_bal_in
        type(reference_t) :: reference

        integer :: opts, ierr
        character(10), parameter :: keys(8) = [character(10) :: 'sys', 'qmc', 'dmqmc', 'ipdmqmc', 'operators', 'rdm', &
                                                                'restart', 'reference']

        lua_state = flu_copyptr(L)
        call get_sys_t(lua_state, sys)

        opts = aot_table_top(lua_state)
        call read_qmc_in(lua_state, opts, qmc_in)
        call read_dmqmc_in(lua_state, sys%basis%nbasis, opts, sys%system, dmqmc_in, dmqmc_estimates%subsys_info)

        ! We are required to handle the tensor length ourselves.
        sys%basis%tensor_label_len = 2*sys%basis%string_len

        if (doing_dmqmc_calc(dmqmc_energy_squared)) then
            ! Create info no longer set in init_real_space.
            allocate(sys%real_lattice%next_nearest_orbs(sys%basis%nbasis,sys%basis%nbasis), stat=ierr)
            call check_allocate('sys%real_lattice%next_nearest_orbs',sys%basis%nbasis*sys%basis%nbasis,ierr)
            call create_next_nearest_orbs(sys%basis, sys%real_lattice)
        end if

        ! [todo] - do spin polarisation in system setup.
        ! I'm sure we can handle this more gracefully than hack it in here...!
        if (dmqmc_in%all_spin_sectors) then
            if (sys%system /= read_in .and. sys%system /= ueg) sys%Ms = sys%lattice%nsites
            call set_spin_polarisation(sys%basis%nbasis, sys)
            if (sys%system == heisenberg) sys%max_number_excitations = sys%lattice%nsites/2
        else
            call set_spin_polarisation(sys%basis%nbasis, sys)
        end if

        ! Now system initialisation is complete (boo), act on the other options.
        call read_restart_in(lua_state, opts, restart_in)
        call read_reference_t(lua_state, opts, sys, reference)
        ! load balancing not implemented in DMQMC--just use defaults.
        call warn_unused_args(lua_state, keys, opts)
        call aot_table_close(lua_state, opts)

        calc_type = dmqmc_calc
        call init_proc_pointers(sys, qmc_in, reference, dmqmc_in)
        call do_dmqmc(sys, qmc_in, dmqmc_in, dmqmc_estimates, restart_in, load_bal_in, reference)

        sys%basis%tensor_label_len = sys%basis%string_len

        nresult = 0

    end function lua_dmqmc

    ! --- table-derived type wrappers ---

    subroutine read_fci_in(lua_state, opts, basis, fci_in, use_sparse_hamil)

        ! Read in the fci and (optionally) lanczos tables to a fci_in_t object.

        ! fci = {
        !     write_hamiltonian = true/false,
        !     hamiltonian_file = filename,
        !     write_determinants = true/false,
        !     determinant_file = filename,
        !     write_nwfns = M,
        !     wfn_file = filename,
        !     nanalyse = N,
        !     blacs_block_size = block_size,
        !     rdm = { ... }, -- L-d vector containing the sites to include in subsystem A.
        ! }
        ! lanczos = {
        !     neigv = N,
        !     nbasis = M,
        !     direct = true/false,
        !     sparse = true/false, -- default true
        ! }

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !    opts: handle for the table containing the fci and (optionally) the lanczos table(s).
        !    basis: information about the single-particle basis set of the system.
        ! Out:
        !    fci_in: fci_in_t object containing generic fci/lanczos input options.
        !    use_sparse_hamil: should the Hamiltonian be stored in a sparse format?

        use flu_binding, only: flu_State
        use aot_table_module, only: aot_get_val, aot_exists, aot_table_open, aot_table_close
        use aot_vector_module, only: aot_get_val

        use lua_hande_utils, only: warn_unused_args
        use basis_types, only: basis_t
        use checking, only: check_allocate
        use fci_utils, only: fci_in_t

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        type(basis_t), intent(in) :: basis
        type(fci_in_t), intent(inout) :: fci_in
        logical, intent(out) :: use_sparse_hamil

        integer :: fci_table, err, fci_nrdms
        integer, allocatable :: err_arr(:)

        character(18), parameter :: fci_keys(9) = [character(18) :: 'write_hamiltonian', 'hamiltonian_file', &
                                                                    'write_determinants', 'determinant_file', 'write_wfns', &
                                                                    'wfn_file', 'nanalyse', 'blacs_block_size', 'rdm']
        character(6), parameter :: lanczos_keys(4) = [character(6) :: 'neigv', 'nbasis', 'direct', 'sparse']

        if (aot_exists(lua_state, opts, 'fci')) then
            call aot_table_open(lua_state, opts, fci_table, 'fci')

            ! Optional arguments
            call aot_get_val(fci_in%write_hamiltonian, err, lua_state, fci_table, 'write_hamiltonian')
            call aot_get_val(fci_in%hamiltonian_file, err, lua_state, fci_table, 'hamiltonian_file')
            call aot_get_val(fci_in%write_determinants, err, lua_state, fci_table, 'write_determinants')
            call aot_get_val(fci_in%determinant_file, err, lua_state, fci_table, 'determinant_file')
            call aot_get_val(fci_in%print_fci_wfn, err, lua_state, fci_table, 'write_nwfns')
            call aot_get_val(fci_in%print_fci_wfn_file, err, lua_state, fci_table, 'wfn_file')
            call aot_get_val(fci_in%analyse_fci_wfn, err, lua_state, fci_table, 'nanalyse')
            call aot_get_val(fci_in%block_size, err, lua_state, fci_table, 'blacs_block_size')

            ! Optional arguments requiring special care.
            if (aot_exists(lua_state, fci_table, 'rdm')) then
                ! Currently restricted to one RDM in a single FCI calculation.
                fci_nrdms = 1
                allocate(fci_in%subsys_info(fci_nrdms), stat=err)
                call check_allocate('fci_in%subsys_info', fci_nrdms, err)
                call aot_get_val(fci_in%subsys_info(fci_nrdms)%subsystem_A, err_arr, basis%nbasis, lua_state, fci_table, 'rdm')
                fci_in%subsys_info(fci_nrdms)%A_nsites = size(fci_in%subsys_info(fci_nrdms)%subsystem_A)
            end if

            call warn_unused_args(lua_state, fci_keys, fci_table)

            call aot_table_close(lua_state, fci_table)

        end if

        ! Lanczos table: optional and indicates doing a Lanczos calculation.
        if (aot_exists(lua_state, opts, 'lanczos')) then
            call aot_table_open(lua_state, opts, fci_table, 'lanczos')
            call aot_get_val(fci_in%nlanczos_eigv, err, lua_state, fci_table, 'neigv')
            call aot_get_val(fci_in%lanczos_string_len, err, lua_state, fci_table, 'nbasis')
            call aot_get_val(fci_in%direct_lanczos, err, lua_state, fci_table, 'direct')
            call aot_get_val(use_sparse_hamil, err, lua_state, fci_table, 'sparse', default=.true.)
            call warn_unused_args(lua_state, lanczos_keys, fci_table)
            call aot_table_close(lua_state, fci_table)
        end if

    end subroutine read_fci_in

    subroutine read_hilbert_args(lua_state, opts, nel, nattempts, ncycles, ex_level, ref_det, rng_seed)

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !    opts: handle to the table which is input to the Lua hilbert_space
        !        routine.
        !    nel: The number of electrons in the system.
        ! Out:
        !    nattempts: number of samples, i.e. number of random determinants, to use each cycle.
        !    ncycles: number of Monte Carlo cycles to perform.
        !    ex_level: maximum excitation level relative to the reference
        !        determinant to include in the Hilbert space.
        !    ref_det: reference determinant.  If not supplied by the user then
        !        this will be deallocated on output.
        !    rng_seed: seed to initialise the random number generator.

        use flu_binding, only: flu_State
        use aot_table_module, only: aot_get_val, aot_exists, aot_table_open, aot_table_close
        use aot_vector_module, only: aot_get_val

        use parallel, only: parent
        use errors, only: stop_all
        use lua_hande_utils, only: warn_unused_args, get_rng_seed

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts, nel
        integer, intent(out) :: nattempts, ncycles, ex_level, rng_seed
        integer, allocatable :: ref_det(:)

        integer :: hilbert_table, err
        integer, allocatable :: err_arr(:)
        character(10), parameter :: keys(6) = [character(10) :: 'sys', 'nattempts', 'ncycles', 'ex_level', 'reference', 'rng_seed']

        if (.not. aot_exists(lua_state, opts, 'hilbert') .and. parent) &
            call stop_all('read_hilbert_args','"hilbert" table not present.')
        call aot_table_open(lua_state, opts, hilbert_table, 'hilbert')

        call aot_get_val(nattempts, err, lua_state, hilbert_table, 'nattempts')
        if (err /= 0 .and. parent) call stop_all('read_hilbert_args', 'Number of attempts not supplied.')
        call aot_get_val(ncycles, err, lua_state, hilbert_table, 'ncycles', default=20)
        call aot_get_val(ex_level, err, lua_state, hilbert_table, 'ex_level', default=-1)
        call aot_get_val(ref_det, err_arr, nel, lua_state, hilbert_table, key='reference')
        call get_rng_seed(lua_state, hilbert_table, rng_seed)

        call warn_unused_args(lua_state, keys, hilbert_table)
        call aot_table_close(lua_state, hilbert_table)

    end subroutine read_hilbert_args

    subroutine read_canonical_energy_args(lua_state, opts, fermi_temperature, beta, nattempts, ncycles, rng_seed, all_spin_sectors)

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !    opts: handle to the table which is input to the Lua canonical_energy_estimate
        !        routine.
        ! Out:
        !    fermi_temperature: if true, rescale beta as the inverse reduced
        !        temperature.
        !    beta: target temperature.
        !    nattempts: number of samples to use each cycle.
        !    ncycles: number of Monte Carlo cycles to perform.
        !    rng_seed: seed to initialise the random number generator.
        !    all_spin_sectors: True if averaging over spin.

        use flu_binding, only: flu_State
        use aot_table_module, only: aot_get_val, aot_exists, aot_table_open, aot_table_close

        use const, only: p
        use parallel, only: parent
        use errors, only: stop_all
        use lua_hande_utils, only: warn_unused_args, get_rng_seed

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        logical, intent(out) :: fermi_temperature, all_spin_sectors
        integer, intent(out) :: nattempts, ncycles, rng_seed
        real(p), intent(out) :: beta

        integer :: canonical_energy_table, err
        character(17), parameter :: keys(6) = [character(17) :: 'nattempts', 'ncycles', 'beta', &
                                                                'fermi_temperature', 'rng_seed', 'all_spin_sectors']

        if (.not. aot_exists(lua_state, opts, 'canonical_energy') .and. parent) &
            call stop_all('read_canonical_energy_args','"canonical_energy" table not present.')
        call aot_table_open(lua_state, opts, canonical_energy_table, 'canonical_energy')

        call aot_get_val(nattempts, err, lua_state, canonical_energy_table, 'nattempts')
        if (err /= 0 .and. parent) call stop_all('read_canonical_energy_args', 'nattempts: Number of attempts/cycle not supplied.')
        call aot_get_val(ncycles, err, lua_state, canonical_energy_table, 'ncycles')
        if (err /= 0 .and. parent) call stop_all('read_canonical_energy_args', 'ncycles: Number of cycles not supplied.')

        call aot_get_val(beta, err, lua_state, canonical_energy_table, 'beta')
        if (err /= 0 .and. parent) call stop_all('read_canonical_energy_args', 'beta: target temperature not supplied.')
        call aot_get_val(fermi_temperature, err, lua_state, canonical_energy_table, 'fermi_temperature', default=.false.)

        call get_rng_seed(lua_state, canonical_energy_table, rng_seed)
        call aot_get_val(all_spin_sectors, err, lua_state, canonical_energy_table, 'all_spin_sectors', default=.false.)

        call warn_unused_args(lua_state, keys, canonical_energy_table)
        call aot_table_close(lua_state, canonical_energy_table)

    end subroutine read_canonical_energy_args

    subroutine read_qmc_in(lua_state, opts, qmc_in, short)

        ! Read in a qmc table to a qmc_in_t object.

        ! qmc = {
        !     tau = tau,                                  -- required
        !     init_pop = N,                               -- required
        !     mc_cycles = ncycles,                        -- required
        !     nreports = nreport,                         -- required
        !     state_size = walker_length,                 -- required
        !     spawned_state_size = spawned_walker_length, -- required
        !     rng_seed = seed,
        !     target_population = pop,
        !     real_amplitudes = true/false,
        !     real_amplitude_force_32 = true/false,
        !     spawn_cutoff = cutoff,
        !     no_renorm = true/false,
        !     tau_search = true/false,
        !     pattempt_single = prob,
        !     pattempt_double = prob,
        !     initial_shift = shift,
        !     shift_damping = damp_factor,
        !     initiator = true/false,
        !     initiator_threshold = pop,
        !     use_mpi_barriers = true/false,
        !     vary_shift_from = shift or "proje",
        ! }

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !    opts: handle for the table containing the qmc table.
        !    short (optional): if true, then don't read in otherwise critical
        !       options not required for simple_fciqmc.  Default: false.
        ! Out:
        !    qmc_in: qmc_in_t object containing generic QMC input options.

        use flu_binding, only: flu_State
        use aot_table_module, only: aot_get_val, aot_exists, aot_table_open, aot_table_close

        use qmc_data, only: qmc_in_t, excit_gen_renorm, excit_gen_no_renorm
        use lua_hande_utils, only: warn_unused_args, get_rng_seed
        use parallel, only: parent
        use errors, only: stop_all, warning

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        type(qmc_in_t), intent(out) :: qmc_in
        logical, intent(in), optional :: short

        integer :: qmc_table, err
        character(len=10) :: str
        logical :: skip, no_renorm

        character(23), parameter :: keys(22) = [character(23) :: 'tau', 'init_pop', 'mc_cycles', 'nreports', 'state_size', &
                                                                 'spawned_state_size', 'rng_seed', 'target_population', &
                                                                 'real_amplitudes', 'spawn_cutoff', 'no_renorm', 'tau_search', &
                                                                 'real_amplitude_force_32', &
                                                                 'pattempt_single', 'pattempt_double', 'initial_shift', &
                                                                 'shift_damping', 'initiator', 'initiator_threshold', &
                                                                 'use_mpi_barriers', 'vary_shift_from', 'excit_gen']

        if (present(short)) then
            skip = short
        else
            skip = .false.
        end if

        if (.not. aot_exists(lua_state, opts, 'qmc') .and. parent) call stop_all('read_qmc_in','"qmc" table not present.')

        call aot_table_open(lua_state, opts, qmc_table, 'qmc')

        ! Required arguments
        call aot_get_val(qmc_in%tau, err, lua_state, qmc_table, 'tau')
        if (err /= 0 .and. parent) call stop_all('read_qmc_in', 'tau not set.')
        call aot_get_val(qmc_in%D0_population, err, lua_state, qmc_table, 'init_pop')
        if (err /= 0 .and. parent) call stop_all('read_qmc_in', 'init_pop not set.')
        call aot_get_val(qmc_in%ncycles, err, lua_state, qmc_table, 'mc_cycles')
        if (err /= 0 .and. parent) call stop_all('read_qmc_in', 'mc_cycles not set.')
        call aot_get_val(qmc_in%nreport, err, lua_state, qmc_table, 'nreports')
        if (err /= 0 .and. parent) call stop_all('read_qmc_in', 'nreports not set.')
        if (.not. skip) then
            call aot_get_val(qmc_in%walker_length, err, lua_state, qmc_table, 'state_size')
            if (err /= 0 .and. parent) call stop_all('read_qmc_in', 'state_size not set.')
            call aot_get_val(qmc_in%spawned_walker_length, err, lua_state, qmc_table, 'spawned_state_size')
            if (err /= 0 .and. parent) call stop_all('read_qmc_in', 'spawned_state_size not set.')
        end if

        ! Optional arguments (defaults set in derived type).
        call aot_get_val(qmc_in%real_amplitudes, err, lua_state, qmc_table, 'real_amplitudes')
        call aot_get_val(qmc_in%real_amplitude_force_32, err, lua_state, qmc_table, 'real_amplitude_force_32')
        call aot_get_val(qmc_in%spawn_cutoff, err, lua_state, qmc_table, 'spawn_cutoff')
        call aot_get_val(qmc_in%pattempt_single, err, lua_state, qmc_table, 'pattempt_single')
        call aot_get_val(qmc_in%pattempt_double, err, lua_state, qmc_table, 'pattempt_double')
        call aot_get_val(qmc_in%tau_search, err, lua_state, qmc_table, 'tau_search')
        call aot_get_val(qmc_in%initial_shift, err, lua_state, qmc_table, 'initial_shift')
        call aot_get_val(qmc_in%shift_damping, err, lua_state, qmc_table, 'shift_damping')
        call aot_get_val(qmc_in%target_particles, err, lua_state, qmc_table, 'target_population')
        call aot_get_val(qmc_in%initiator_approx, err, lua_state, qmc_table, 'initiator')
        call aot_get_val(qmc_in%initiator_pop, err, lua_state, qmc_table, 'initiator_threshold')
        call aot_get_val(qmc_in%use_mpi_barriers, err, lua_state, qmc_table, 'use_mpi_barriers')

        if (aot_exists(lua_state, qmc_table, 'no_renorm')) then
            call aot_get_val(no_renorm, err, lua_state, qmc_table, 'no_renorm')
            call warning('read_qmc_in', 'no_renorm is deprecated.  Please use excit_gen instead.')
            if (no_renorm) then
                qmc_in%excit_gen = excit_gen_no_renorm
            else
                qmc_in%excit_gen = excit_gen_renorm
            end if
        end if

        if (aot_exists(lua_state, qmc_table, 'excit_gen')) then
            call aot_get_val(str, err, lua_state, qmc_table, 'excit_gen')
            select case(str)
            case('renorm')
                qmc_in%excit_gen = excit_gen_renorm
            case('no_renorm')
                qmc_in%excit_gen = excit_gen_no_renorm
            case default
                call stop_all('read_qmc_in', 'Invalid excit_gen setting: '//trim(str))
            end select
        end if

        ! If user sets initial shift and vary_shift_from, assume they know what
        ! they're doing.  Otherwise, vary the shift from the initial shift
        ! value.
        qmc_in%vary_shift_from = qmc_in%initial_shift

        ! Optional arguments requiring special care.
        call get_rng_seed(lua_state, qmc_table, qmc_in%seed)
        if (aot_exists(lua_state, qmc_table, 'vary_shift_from')) then
            call aot_get_val(qmc_in%vary_shift_from, err, lua_state, qmc_table, 'vary_shift_from')
            if (err /= 0) then
                ! Perhaps passed in a string?
                call aot_get_val(str, err, lua_state, qmc_table, 'vary_shift_from')
                qmc_in% vary_shift_from_proje = trim(str) == 'proje'
            end if
        end if

        call warn_unused_args(lua_state, keys, qmc_table)

        call aot_table_close(lua_state, qmc_table)

    end subroutine read_qmc_in

    subroutine read_fciqmc_in(lua_state, opts, fciqmc_in)

        ! Read in an fciqmc table (if it exists) to an fciqmc_in object.

        ! fciqmc = {
        !     non_blocking_comm = true/false,
        !     load_balancing = true/false,
        !     init_spin_inverse_reference_det = true/false,
        !     guiding_function = 'none'/'neel_singlet',
        !     trial_function = 'single_basis'/'neel_singlet',
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

        use qmc_data, only: fciqmc_in_t, neel_singlet, neel_singlet_guiding
        use lua_hande_utils, only: warn_unused_args
        use parallel, only: parent
        use errors, only: stop_all

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        type(fciqmc_in_t), intent(out) :: fciqmc_in

        integer :: fciqmc_table, ref_det, err
        character(len=12) :: str
        logical :: ref_det_flag
        character(31), parameter :: keys(6) = [character(31) :: 'non_blocking_comm', 'load_balancing', 'guiding_function', &
                                                                'init_spin_inverse_reference_det', 'trial_function', &
                                                                'select_reference_det']

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
                    call aot_get_val(ref_det_flag, err, lua_state, fciqmc_table, 'select_reference_det', default=.false.)
                    if (ref_det_flag) fciqmc_in%select_ref_det_every_nreports = 20
                else
                    call aot_get_val(fciqmc_in%select_ref_det_every_nreports, err, lua_state, ref_det, 'update_every', default=20)
                    call aot_get_val(fciqmc_in%ref_det_factor, err, lua_state, ref_det, 'pop_factor')
                    call aot_table_close(lua_state, ref_det)
                end if
            end if
            if (aot_exists(lua_state, fciqmc_table, 'guiding_function')) then
                call aot_get_val(str, err, lua_state, fciqmc_table, 'guiding_function')
                select case(trim(str))
                case('none')
                    ! use default
                case('neel_singlet')
                    fciqmc_in%guiding_function = neel_singlet_guiding
                    fciqmc_in%trial_function = neel_singlet
                case default
                    if (parent) call stop_all('read_fciqmc_in', 'Unknown guiding function')
                end select
            end if
            if (aot_exists(lua_state, fciqmc_table, 'trial_function')) then
                call aot_get_val(str, err, lua_state, fciqmc_table, 'trial_function')
                select case(trim(str))
                case('none')
                    ! use default
                case('neel_singlet')
                    fciqmc_in%trial_function = neel_singlet
                case default
                    if (parent) call stop_all('read_fciqmc_in', 'Unknown trial wavefunction')
                end select
            end if

            call warn_unused_args(lua_state, keys, fciqmc_table)

            call aot_table_close(lua_state, fciqmc_table)

        end if

    end subroutine read_fciqmc_in

    subroutine read_semi_stoch_in(lua_state, opts, qmc_in, semi_stoch_in)

        ! Read in an semi_stoch table (if it exists) to an semi_stoch_in object.

        ! semi_stoch = {
        !     space = 'read'/'high',            -- required if semi_stoch table exists.
        !     size = target_pop,                -- required if space == 'high'
        !     start_iteration = mc_cycles,
        !     shift_start_iteration = mc_cycle, -- overrides start_iteration
        !     separate_annihilation = true/false,
        !     write = true/false/id,
        !     read = id,
        ! }

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        !    qmc_in: qmc_in_t object containing qmc input options.  If semi-stochastic
        !        options are present, then qmc_in%real_amplitudes is forced to be true.
        ! In:
        !    opts: handle for the table containing the semi_stoch table.
        ! Out:
        !    semi_stoch_in: semi_stoch_in_t object containing semi-stochastic-specific input options.

        use flu_binding, only: flu_State
        use aot_table_module, only: aot_get_val, aot_exists, aot_table_open, aot_table_close
        use parallel, only: parent
        use errors, only: stop_all, warning

        use qmc_data, only: qmc_in_t, semi_stoch_in_t, high_pop_determ_space, read_determ_space, &
                            semi_stoch_combined_annihilation, semi_stoch_separate_annihilation
        use lua_hande_utils, only: warn_unused_args, get_flag_and_id

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        type(qmc_in_t), intent(inout) :: qmc_in
        type(semi_stoch_in_t), intent(out) :: semi_stoch_in

        integer :: semi_stoch_table, err
        character(len=10) :: str
        logical :: separate_annihilation
        character(21), parameter :: keys(8) = [character(21) :: 'space', 'size', 'start_iteration', 'separate_annihilation', &
                                                                'shift_start_iteration', 'write_determ_space', 'write', 'read']

        if (aot_exists(lua_state, opts, 'semi_stoch')) then

            qmc_in%real_amplitudes = .true.

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
                    if (err /= 0 .and. parent) call stop_all('read_semi_stoch_in', 'Target space size not given.')
                case default
                    if (parent) call stop_all('read_semi_stoch_in', 'Unknown semi stochastic space type')
                end select
            else
                if (parent) call stop_all('read_semi_stoch_in', 'Deterministic space not specified.')
            end if

            ! Optional arguments (defaults set in derived type).
            call aot_get_val(semi_stoch_in%start_iter, err, lua_state, semi_stoch_table, 'start_iteration')
            call aot_get_val(separate_annihilation, err, lua_state, semi_stoch_table, 'separate_annihilation', default=.true.)
            call aot_get_val(semi_stoch_in%read_id, err, lua_state, semi_stoch_table, 'read')
            call get_flag_and_id(lua_state, semi_stoch_table, semi_stoch_in%write_determ_space, semi_stoch_in%write_id, 'write')
            if (separate_annihilation) then
                semi_stoch_in%projection_mode = semi_stoch_separate_annihilation
            else
                semi_stoch_in%projection_mode = semi_stoch_combined_annihilation
            end if

            ! Optional arguments requiring special care.
            if (aot_exists(lua_state, semi_stoch_table, 'write_determ_space')) then
                call aot_get_val(semi_stoch_in%write_determ_space, err, lua_state, semi_stoch_table, 'write_determ_space')
                call warning('read_semi_stoch_in', 'write_determ_space is deprecated.  Please use write instead.')
            end if
            if (aot_exists(lua_state, semi_stoch_table, 'shift_start_iteration')) then
                semi_stoch_in%start_iter = huge(0) ! fixed up once the shift comes on...
                call aot_get_val(semi_stoch_in%shift_iter, err, lua_state, semi_stoch_table, 'shift_start_iteration')
            end if

            call warn_unused_args(lua_state, keys, semi_stoch_table)

            call aot_table_close(lua_state, semi_stoch_table)

        end if

    end subroutine read_semi_stoch_in

    subroutine read_ccmc_in(lua_state, opts, ccmc_in)

        ! Read in an ccmc table (if it exists) to an ccmc_in object.

        ! ccmc = {
        !     move_frequency = frequency,
        !     cluster_multispawn_threshold = threshold,
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

        integer :: ccmc_table, err
        character(28), parameter :: keys(4) = [character(28) :: 'move_frequency', 'cluster_multispawn_threshold', &
                                                                'full_non_composite', 'linked']

        if (aot_exists(lua_state, opts, 'ccmc')) then

            call aot_table_open(lua_state, opts, ccmc_table, 'ccmc')

            ! Optional arguments (defaults set in derived type).
            call aot_get_val(ccmc_in%move_freq, err, lua_state, ccmc_table, 'move_frequency')
            call aot_get_val(ccmc_in%cluster_multispawn_threshold, err, lua_state, ccmc_table, 'cluster_multispawn_threshold')
            call aot_get_val(ccmc_in%full_nc, err, lua_state, ccmc_table, 'full_non_composite')
            call aot_get_val(ccmc_in%linked, err, lua_state, ccmc_table, 'linked')

            call warn_unused_args(lua_state, keys, ccmc_table)

            call aot_table_close(lua_state, ccmc_table)

        end if

    end subroutine read_ccmc_in

    subroutine read_dmqmc_in(lua_state, nbasis, opts, system_name, dmqmc_in, subsys_info)

        ! Read in DMQMC-related input options from various tables to a dmqmc_in_t object.

        ! dmqmc = {
        !     replica_tricks = true/false,
        !     fermi_temperature = true/false,
        !     all_sym_sectors = true/false,
        !     all_spin_sectors = true/false,
        !     beta_loops = Nb,
        !     sampling_weights = { ... },
        !     vary_weights = N,
        !     find_weights = true/false,
        !     symmetrize = true/false,
        !     initiator_level = ilevel,
        ! }
        ! ipdmqmc = { -- sets propagate_to_beta to true
        !     initial_beta = b,
        !     initial_matrix = 'free_electron'/'hartree_fock',
        !     grand_canonical_initialisation = true/false,
        !     metropolis_attempts = nattempts,
        !     symmetric = true/false,
        ! }
        ! operators = {
        !     renyi2 = true/false,
        !     energy = true/false,
        !     energy2 = true/false,
        !     staggered_magnetisation = true/false,
        !     correlation = { ... }
        !     excit_dist = true/false,
        !     excit_dist_start = iteration,
        !     kinetic_energy = true/false,
        !     H0_energy = true/false,
        !     potential_energy = true/false,
        !     HI_energy = true/false,
        ! }
        ! rdm = {
        !     spawned_state_size = X,
        !     rdms = { { ... } , { ... } , { ... } ... }
        !     ground_state = true/false,
        !     ground_state_start = iteration,
        !     instantaneous = true/false,
        !     write = true/false,
        !     concurrence = true/false,
        !     von_neumann = true/false,
        !     renyi2 = true/false,
        ! }

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !    opts: handle for the table containing the above tables.
        !    nbasis: number of single-particle basis functions in the basis set of the system.
        !    system_name: system being studied.
        ! Out:
        !    dmqmc_in: dmqmc_in_t object containing DMQMC input options.
        !    subsys_info: array of subsys_t objects containing the subsystem
        !       for each RDM of interest. Unallocated if no RDMs are specified.

        use flu_binding, only: flu_State
        use aot_table_module, only: aot_exists, aot_table_open, aot_table_close, aot_get_val, aot_table_length
        use aot_vector_module, only: aot_get_val

        use calc
        use dmqmc_data, only: dmqmc_in_t, subsys_t, hartree_fock_dm, free_electron_dm
        use checking, only: check_allocate
        use parallel, only: parent
        use errors, only: stop_all
        use lua_hande_utils, only: warn_unused_args
        use system, only: heisenberg

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: nbasis, opts, system_name
        type(dmqmc_in_t), intent(out) :: dmqmc_in
        type(subsys_t), intent(out), allocatable :: subsys_info(:)

        integer :: table, subtable, err, i
        integer, allocatable :: err_arr(:)
        logical :: op
        character(len=13) :: str

        character(30), parameter :: dmqmc_keys(11) = [character(30) :: 'replica_tricks', 'fermi_temperature', 'all_sym_sectors', &
                                                                      'all_spin_sectors', 'beta_loops', 'sampling_weights',      &
                                                                      'find_weights', 'find_weights_start', 'symmetrize',        &
                                                                      'vary_weights', 'initiator_level']
        character(30), parameter :: ip_keys(5)    = [character(30) :: 'initial_beta', 'initial_matrix',                          &
                                                                      'grand_canonical_initialisation', 'metropolis_attempts',   &
                                                                      'symmetric']
        character(30), parameter :: op_keys(10)    = [character(30) :: 'renyi2', 'energy', 'energy2', 'staggered_magnetisation',  &
                                                                       'correlation', 'excit_dist', 'kinetic_energy',             &
                                                                       'H0_energy', 'potential_energy', 'HI_energy']
        character(30), parameter :: rdm_keys(9)   = [character(30) :: 'spawned_state_size', 'rdms', 'ground_state',              &
                                                                      'ground_rdm_start', 'instantaneous', 'write',              &
                                                                      'concurrence', 'von_neumann', 'renyi2']

        dmqmc_calc_type = 0

        ! All optional and straightfoward except the vector quantities.

        if (aot_exists(lua_state, opts, 'dmqmc')) then
            call aot_table_open(lua_state, opts, table, 'dmqmc')
            call aot_get_val(dmqmc_in%replica_tricks, err, lua_state, table, 'replica_tricks')
            call aot_get_val(dmqmc_in%fermi_temperature, err, lua_state, table, 'fermi_temperature')
            call aot_get_val(dmqmc_in%all_sym_sectors, err, lua_state, table, 'all_sym_sectors')
            call aot_get_val(dmqmc_in%all_spin_sectors, err, lua_state, table, 'all_spin_sectors')
            call aot_get_val(dmqmc_in%beta_loops, err, lua_state, table, 'beta_loops')
            call aot_get_val(dmqmc_in%find_weights, err, lua_state, table, 'find_weights')
            call aot_get_val(dmqmc_in%find_weights_start, err, lua_state, table, 'find_weights_start')
            call aot_get_val(dmqmc_in%half_density_matrix, err, lua_state, table, 'symmetrize')
            if (aot_exists(lua_state, table, 'sampling_weights')) then
                dmqmc_in%weighted_sampling = .true.
                ! Certainly can't have more excitation levels than basis functions, so that's a handy upper-limit.
                call aot_get_val(dmqmc_in%sampling_probs, err_arr, nbasis, lua_state, table, 'sampling_weights')
            end if
            if (aot_exists(lua_state, table, 'initiator_level')) then
                call aot_get_val(dmqmc_in%initiator_level, err, lua_state, table, 'initiator_level')
            end if
            dmqmc_in%vary_weights = aot_exists(lua_state, table, 'vary_weights')
            call aot_get_val(dmqmc_in%finish_varying_weights, err, lua_state, table, 'vary_weights')
            call warn_unused_args(lua_state, dmqmc_keys, table)
            call aot_table_close(lua_state, table)
        end if

        ! If using the 'find_weights' option then weighted sampling should
        ! always be used, even if not asked for via a separate input.
        if (dmqmc_in%find_weights) dmqmc_in%weighted_sampling = .true.

        if (aot_exists(lua_state, opts, 'ipdmqmc')) then
            dmqmc_in%propagate_to_beta = .true.
            call aot_table_open(lua_state, opts, table, 'ipdmqmc')
            call aot_get_val(dmqmc_in%init_beta, err, lua_state, table, 'initial_beta')
            if (aot_exists(lua_state, table, 'initial_matrix')) then
                call aot_get_val(str, err, lua_state, table, 'initial_matrix')
                select case(trim(str))
                case('free_electron')
                    dmqmc_in%initial_matrix = free_electron_dm
                case('hartree_fock')
                    dmqmc_in%initial_matrix = hartree_fock_dm
                case default
                    if (parent) call stop_all('read_dmqmc_in', 'Unknown  inital density matrix')
                end select
            end if
            call aot_get_val(dmqmc_in%symmetric, err, lua_state, table, 'symmetric', default=.true.)
            call aot_get_val(dmqmc_in%grand_canonical_initialisation, err, lua_state, table, 'grand_canonical_initialisation')
            call aot_get_val(dmqmc_in%metropolis_attempts, err, lua_state, table, 'metropolis_attempts')
            call warn_unused_args(lua_state, ip_keys, table)
            call aot_table_close(lua_state, table)
        end if

        if (aot_exists(lua_state, opts, 'operators')) then
            call aot_table_open(lua_state, opts, table, 'operators')
            call aot_get_val(op, err, lua_state, table, 'renyi2', default=.false.)
            if (op) dmqmc_calc_type = dmqmc_calc_type + dmqmc_full_r2
            call aot_get_val(op, err, lua_state, table, 'energy', default=.false.)
            if (op) dmqmc_calc_type = dmqmc_calc_type + dmqmc_energy
            call aot_get_val(op, err, lua_state, table, 'energy2', default=.false.)
            if (op) dmqmc_calc_type = dmqmc_calc_type + dmqmc_energy_squared
            call aot_get_val(op, err, lua_state, table, 'staggered_magnetisation', default=.false.)
            if (op) dmqmc_calc_type = dmqmc_calc_type + dmqmc_staggered_magnetisation
            call aot_get_val(op, err, lua_state, table, 'kinetic_energy', default=.false.)
            if (op) dmqmc_calc_type = dmqmc_calc_type + dmqmc_kinetic_energy
            call aot_get_val(op, err, lua_state, table, 'H0_energy', default=.false.)
            if (op) dmqmc_calc_type = dmqmc_calc_type + dmqmc_H0_energy
            call aot_get_val(op, err, lua_state, table, 'potential_energy', default=.false.)
            if (op) dmqmc_calc_type = dmqmc_calc_type + dmqmc_potential_energy
            call aot_get_val(op, err, lua_state, table, 'HI_energy', default=.false.)
            if (op) dmqmc_calc_type = dmqmc_calc_type + dmqmc_HI_energy
            if (aot_exists(lua_state, table, 'correlation')) then
                dmqmc_calc_type = dmqmc_calc_type + dmqmc_correlation
                call aot_get_val(dmqmc_in%correlation_sites, err_arr, nbasis, lua_state, table, 'correlation')
            end if
            call aot_get_val(dmqmc_in%calc_excit_dist, err, lua_state, table, 'excit_dist')
            call warn_unused_args(lua_state, op_keys, table)
            call aot_table_close(lua_state, table)
        end if

        if (aot_exists(lua_state, opts, 'rdm')) then
            call aot_table_open(lua_state, opts, table, 'rdm')

            ! instantaneous is optional but if present, spawned_state_size must
            ! be.
            call aot_get_val(dmqmc_in%rdm%calc_inst_rdm, err, lua_state, table, 'instantaneous')

            ! If doing rdms, must specify them and the spawned state size.
            call aot_get_val(dmqmc_in%rdm%spawned_length, err, lua_state, table, 'spawned_state_size')
            if (dmqmc_in%rdm%calc_inst_rdm .and. err /= 0) then
                if (parent) call stop_all('read_dmqmc_in', 'spawned_state_size not specified in rdm table.')
            end if
            if (aot_exists(lua_state, table, 'rdms')) then
                dmqmc_in%rdm%doing_rdm = .true.
                call aot_table_open(lua_state, table, subtable, 'rdms')
                dmqmc_in%rdm%nrdms = aot_table_length(lua_state, subtable)
                allocate(subsys_info(dmqmc_in%rdm%nrdms), stat=err)
                call check_allocate('subsys_info', dmqmc_in%rdm%nrdms, err)
                do i = 1, dmqmc_in%rdm%nrdms
                    call aot_get_val(subsys_info(i)%subsystem_A, err_arr, nbasis, lua_state, subtable, pos=i)
                    subsys_info(i)%A_nsites = size(subsys_info(i)%subsystem_A)
                end do
            else
                if (parent) call stop_all('read_dmqmc_in', 'rdms not specified in rdm table.')
            end if

            ! Optional arguments.
            call aot_get_val(dmqmc_in%rdm%calc_ground_rdm, err, lua_state, table, 'ground_state')
            call aot_get_val(dmqmc_in%start_av_rdm, err, lua_state, table, 'ground_rdm_start')
            call aot_get_val(dmqmc_in%rdm%output_rdm, err, lua_state, table, 'write')
            call aot_get_val(dmqmc_in%rdm%doing_concurrence, err, lua_state, table, 'concurrence')
            call aot_get_val(dmqmc_in%rdm%doing_vn_entropy, err, lua_state, table, 'von_neumann')
            call aot_get_val(op, err, lua_state, table, 'renyi2', default=.false.)
            if (op) dmqmc_calc_type = dmqmc_calc_type + dmqmc_rdm_r2
            call warn_unused_args(lua_state, rdm_keys, table)
            call aot_table_close(lua_state, table)
        else
            dmqmc_in%rdm%nrdms = 0
        end if
        ! Can't average over spin sectors alone in normal dmqmc calculation.
        if (dmqmc_in%all_spin_sectors .and. .not. dmqmc_in%all_sym_sectors .and. system_name /= heisenberg .and. parent) &
                                            & call stop_all('read_dmqmc_in', 'specified symmetry averaging not imlemented')

    end subroutine read_dmqmc_in

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

        use lua_hande_utils, only: warn_unused_args
        use qmc_data, only: load_bal_in_t

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        type(load_bal_in_t), intent(out) :: load_bal_in

        integer :: load_bal_table, err
        character(12), parameter :: keys(5) = [character(12) :: 'nslots', 'min_pop', 'target', 'max_attempts', 'write']

        if (aot_exists(lua_state, opts, 'load_bal')) then

            call aot_table_open(lua_state, opts, load_bal_table, 'load_bal')

            ! Optional arguments (defaults set in derived type).
            load_bal_in%nslots = 20 ! Default for when load balancing in use
            call aot_get_val(load_bal_in%nslots, err, lua_state, load_bal_table, 'nslots')
            call aot_get_val(load_bal_in%pop, err, lua_state, load_bal_table, 'min_pop')
            call aot_get_val(load_bal_in%percent, err, lua_state, load_bal_table, 'target')
            call aot_get_val(load_bal_in%max_attempts, err, lua_state, load_bal_table, 'max_attempts')
            call aot_get_val(load_bal_in%write_info, err, lua_state, load_bal_table, 'write')

            call warn_unused_args(lua_state, keys, load_bal_table)

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

        use lua_hande_utils, only: warn_unused_args
        use qmc_data, only: reference_t
        use system, only: sys_t

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        type(sys_t), intent(in) :: sys
        type(reference_t), intent(out) :: ref

        integer :: ref_table, err
        integer, allocatable :: err_arr(:)
        character(17), parameter :: keys(3) = [character(17) :: 'det', 'hilbert_space_det', 'ex_level']

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

            call warn_unused_args(lua_state, keys, ref_table)

            call aot_table_close(lua_state, ref_table)

        end if

    end subroutine read_reference_t

    subroutine read_restart_in(lua_state, opts, restart_in)

        ! Read in a restart table (if it exists) to a restart_in_t object.

        ! restart = {
        !    read = true/false/id,
        !    write = true/false/id,
        !    write_shift = true/false/id,
        !    write_frequency = niterations,
        ! }

        ! where id is an integer and sets the id for the restart file of the relevant type.
        ! If an id is provided, then the corresponding flag is also set to true.

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !    opts: handle for the table containing the restart table.
        ! Out:
        !    restart_in: restart_in_t object contianing the input options describing restart files.

        use flu_binding, only: flu_State
        use aot_table_module, only: aot_exists, aot_table_open, aot_table_close, aot_get_val

        use lua_hande_utils, only: warn_unused_args, get_flag_and_id
        use qmc_data, only: restart_in_t

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        type(restart_in_t), intent(out) :: restart_in

        integer :: err, restart_table
        character(15), parameter :: keys(4) = [character(15) :: 'read', 'write', 'write_shift', 'write_frequency']

        if (aot_exists(lua_state, opts, 'restart')) then

            call aot_table_open(lua_state, opts, restart_table, 'restart')

            call aot_get_val(restart_in%write_freq, err, lua_state, restart_table, 'write_frequency')

            associate(r_in=>restart_in)
                call get_flag_and_id(lua_state, restart_table, r_in%read_restart, r_in%read_id, 'read')
                call get_flag_and_id(lua_state, restart_table, r_in%write_restart, r_in%write_id, 'write')
                call get_flag_and_id(lua_state, restart_table, r_in%write_restart_shift, r_in%write_shift_id, 'write_shift')
            end associate

            call warn_unused_args(lua_state, keys, restart_table)

            call aot_table_close(lua_state, restart_table)

        end if

    end subroutine read_restart_in

end module lua_hande_calc
