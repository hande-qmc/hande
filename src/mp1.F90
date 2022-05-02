module mp1

use const, only: p, int_64
use qmc_data, only: excit_gen_renorm

implicit none

! [todo] - document
type mp1_in_t
    real(p) :: D0_norm
    integer :: state_size
    logical :: real_amplitudes = .false.
    real(p) :: spawn_cutoff = 0.01_p
    logical :: even_selection = .false.
end type mp1_in_t



contains

    subroutine sample_mp1_wfn(sys, mp1_in, ref_in, tijab, rng_seed)

        ! Generate the mp1 wavefunction (1+T_2)|HF> exactly and use it as an initial guess for FCIQMC/CCMC
        ! In:
        !   sys: system being studied.
        !   mpi_in: input options relating to mp1.
        !   ref_in: current reference determinant, set in lua_hande_calc.
        !   rng_seed (optional): the integer rng seed used in coarse-graining the wavefunction
        ! Out:
        !   tijab: particle_t object for use in FCIQMC/CCMC calculation

        use const, only: int_p, i0, depsilon

        use checking, only: check_allocate
        use errors, only: stop_all, warning
        use json_out
        use dSFMT_interface, only: dSFMT_t, dSFMT_init, get_rand_close_open
        use parallel

        use system, only: sys_t, copy_sys_spin_info, sys_t_json
        use calc_system_init, only: set_spin_polarisation
        use qmc_data, only: particle_t, reference_t, annihilation_flags_t
        use qmc_data, only: load_bal_state_t, reference_t_json, qmc_in_t, qmc_state_t
        use qmc, only: init_proc_pointers, init_reference, init_excit_gen
        use ccmc_data, only: cluster_t, ex_lvl_dist_t
        use determinants, only: det_info_t, alloc_det_info_t
        use spawn_data, only: spawn_t, proc_map_t, alloc_spawn_t, dealloc_spawn_t

        use annihilation, only: direct_annihilation, remove_unoccupied_dets
        use calc, only: GLOBAL_META, gen_seed
        use load_balancing, only: init_proc_map_t
        use ccmc_utils, only: cumulative_population, find_D0
        use ccmc_selection, only: select_cluster
        use determinants, only: encode_det
        use determinant_decoders, only: decode_det_occ
        use excitations, only: excit_t, create_excited_det, get_excitation_level
        use hamiltonian, only: get_hmatel
        use hamiltonian_data
        use particle_t_utils, only: init_particle_t, dealloc_particle_t
        use proc_pointers, only: decoder_ptr
        use reference_determinant, only: set_reference_det
        use stoch_utils, only: stochastic_round_spawned_particle
        use spawning, only: create_spawned_particle, attempt_to_spawn, assign_particle_processor
        use excit_gens, only: excit_gen_data_t

        type(sys_t), intent(inout) :: sys
        type(mp1_in_t), intent(in) :: mp1_in
        type(reference_t), intent(in), optional :: ref_in
        type(particle_t), intent(out) :: tijab
        integer, intent(in), optional :: rng_seed

        type(spawn_t) :: spawn
        type(reference_t) :: ref
        type(det_info_t) :: cdet
        type(cluster_t) :: cluster
        type(excit_t) :: excit
        type(annihilation_flags_t) :: annihilation_flags
        type(proc_map_t) :: proc_map
        type(dSFMT_t) :: rng
        type(sys_t) :: sys_bak

        integer :: seed, max_nspawned_states, nspaces, state_size, psip_element_size, homo, lumo
        integer :: i, j, a, b, ia, ib, iocc_a, iocc_b
        integer(int_p) :: nspawn
        real(p) :: spawn_cutoff, denom
        type(hmatel_t) :: intgrl, ampl
        real(dp) :: emp2
        integer(i0) :: f(sys%basis%bit_string_len)
        type(json_out_t) :: js
        integer :: io_unit
        type(excit_gen_data_t) :: excit_gen_data
        type(qmc_in_t) :: qmc_in_loc

        io_unit = 6
        if (parent) then
            write (io_unit,'(1X,"Deterministic MP1 wavefunction initialisation")')
            write (io_unit,'(1X, 45("-"))')
        end if

        if (present(rng_seed)) then
            seed = rng_seed
        else
            seed = gen_seed(GLOBAL_META%uuid)
        end if
        call dSFMT_init(seed+iproc, 50000, rng)

        call copy_sys_spin_info(sys, sys_bak)
        call set_spin_polarisation(sys%basis%nbasis, sys)
       
        ! We have to have a qmc_in_t object for init_proc_pointers, as init_reference requires proc_pointers to be properly set. Copied from qmc_data.f90::qmc_in_t
        qmc_in_loc%excit_gen = excit_gen_renorm
        qmc_in_loc%pattempt_single = -1
        qmc_in_loc%pattempt_double = -1
        qmc_in_loc%pattempt_update = .false.
        qmc_in_loc%pattempt_zero_accum_data = .false.
        qmc_in_loc%pattempt_parallel = -1
        call init_proc_pointers(sys, qmc_in_loc, ref, io_unit)
        call init_reference(sys, ref_in, io_unit, ref)
        call init_excit_gen(sys, qmc_in_loc, ref, .false., excit_gen_data)
        
        ! Similarly for particle stores.  We use
        ! * tijab: exact amplitudes on reference and doubles.
        state_size = mp1_in%state_size
        if (state_size < 0) then
            ! Given in MB.  Convert to elements.
            ! bit string and population element.
            psip_element_size = (sys%basis%tensor_label_len*i0_length + int_p_length)/8
            ! data element
#ifdef SINGLE_PRECISION
            psip_element_size = psip_element_size + 4
#else
            psip_element_size = psip_element_size + 8
#endif
            state_size = -int((real(state_size,p)*10**6)/psip_element_size)
        end if

        ! These must be set before passing into init_particle_t
        tijab%nspaces = 1
        tijab%info_size = 0

        call init_particle_t(state_size, 1, sys%basis%tensor_label_len, mp1_in%real_amplitudes, .false., tijab)

        ! annihilation_flags default to 'off'.  Only feature we might be using is real amplitudes.
        annihilation_flags%real_amplitudes = mp1_in%real_amplitudes
        ! Just need the simplest decoder as not doing anything based on the clusters produced beyond setting a bit string...
        decoder_ptr => decode_det_occ

        ! Initialise spawn store.  Use sane defaults for options irrelevant to sampling a wavefunction (i.e. disable initiator, only
        ! one particle space, default processor map, no MPI barriers).
        nspaces = 1
        spawn_cutoff = mp1_in%spawn_cutoff
        call init_proc_map_t(1, proc_map)
        if (.not. mp1_in%real_amplitudes) spawn_cutoff = 0.0_p
        ! Need at least one spawning slot per attempt so give ourselves a generous amount of breathing room for load inbalance.
        max_nspawned_states = 2*ceiling(real(state_size)/nprocs)*nprocs
        call alloc_spawn_t(sys%basis%tensor_label_len, sys%basis%nbasis, nspaces, .false., max_nspawned_states, spawn_cutoff, &
                           tijab%pop_real_factor, proc_map, 7, .false., spawn)
        ! spawn_t now holds all the processor map/load balancing info required so can safely deallocate proc_map.
        deallocate(proc_map%map)

        ! Start with the reference
        tijab%nstates = 1
        tijab%states(:,tijab%nstates) = ref%f0
        tijab%pops(1,tijab%nstates) = nint(mp1_in%D0_norm)*tijab%pop_real_factor
        ! We define the energy zero to be <HF|H|HF>.  The annihilation framework evaluates <D|H|D> for other determinants.
        tijab%dat(1,tijab%nstates) = 0.0_p
        tijab%nparticles(1) = real(sum(abs(tijab%pops(1,:tijab%nstates))),p)/tijab%pop_real_factor

        call alloc_det_info_t(sys, cdet)
        allocate(cluster%excitors(ref%ex_level+2))

        if (parent) then
            call json_object_init(js, tag=.true.)
            call sys_t_json(js, sys)
            call mp1_in_t_json(js, mp1_in)
            call reference_t_json(js, ref, terminal=.true.)
            call json_object_end(js, terminal=.true., tag=.true.)
            write (js%io,'()')
        end if

        if (sys%aufbau_sym) then
            homo = 1
            lumo = sys%basis%nbasis
            associate(bl=>sys%basis%bit_lookup, bfns=>sys%basis%basis_fns)
                do i = 2, sys%basis%nbasis
                    if (btest(ref%f0(bl(2,i)), bl(1,i))) then
                        ! Occupied
                        if (bfns(i)%sp_eigv >= bfns(homo)%sp_eigv) homo = i
                    else
                        ! Unoccupied
                        if (bfns(i)%sp_eigv < bfns(lumo)%sp_eigv) lumo = i
                    end if
                end do
                if (abs(bfns(lumo)%sp_eigv-bfns(homo)%sp_eigv) < depsilon) then
                    call stop_all('sample_mp1_wfn', 'Cannot generate MP1 wavefunction for a system with a zero HOMO-LUMO gap.')
                end if
            end associate
        else
            call warning('sample_mp1_wfn', 'Not using aufbau principle, we do not check that the system is gapless.')
        end if

        ! Generate N_0(1+T2) deterministically.
        excit%nexcit = 2
        
        ! MPI is turned off, the reasons are
        ! 1. The code below tries to find the next excitation by incrementing iocc_a/b, which creates potential
        !   duplicates in excitations on different processors. This leads to different 'exact MP2 energies' with different -np
        ! 2. The code initialises D0_norm number of walkers on the reference on every processor, but we only want 
        !   D0 on one processor in FCIQMC/CCMC
        ! 3. MP2 is not rate limiting

        ! [todo] - OpenMP parallelisation: care is needed to make sure threading and spawning work together (see ccmc.f90)

        if (parent) then
            emp2 = 0.0_dp
            do i = 1, sys%nel
                do j = i+1, sys%nel
                    excit%from_orb = [ref%occ_list0(i), ref%occ_list0(j)]
                    iocc_a = 0
                    iocc_b = 0
                    do ia = 1, sys%basis%nbasis-sys%nel
                        ! Find the next unoccupied orbital (necessary for arbitrary references)
                        ! There should definitely only be (nbasis-nel) virtuals, but we don't assume they're ordered.
                        do
                            ! iocc_a stores the location of the last visited virtual (initialised at 0)
                            a = ia + iocc_a
                            if (.not. btest(ref%f0(sys%basis%bit_lookup(2,a)), sys%basis%bit_lookup(1,a))) exit
                            iocc_a = iocc_a + 1
                        end do
                        do ib = ia, sys%basis%nbasis-sys%nel
                            do
                                b = ib + iocc_b
                                if (.not. btest(ref%f0(sys%basis%bit_lookup(2,b)), sys%basis%bit_lookup(1,b))) exit
                                iocc_b = iocc_b + 1
                            end do

                            excit%to_orb = [a,b]
                            call create_excited_det(sys%basis, ref%f0, excit, f)

                            ! Calculate deterministically the t2 amplitudes and the MP2 energy contributions

                            ! t_{ijab} = <ij||ab> / (e_i + e_j - e_a - e_b).
                            ! NOTE: sp_eigv is not the Hartree-Fock eigenvalues in all cases (e.g. Hubbard model in k-space, UEG!).
                                intgrl = get_hmatel(sys, ref%f0, f)
                                denom = (1.0_p / (sys%basis%basis_fns(i)%sp_eigv + sys%basis%basis_fns(j)%sp_eigv &
                                                - sys%basis%basis_fns(a)%sp_eigv - sys%basis%basis_fns(b)%sp_eigv))
                                ampl = intgrl * denom
                                if (sys%read_in%comp) then
                                    emp2 = emp2 + real(conjg(intgrl%c)*ampl%c)
                                else
                                    emp2 = emp2 + intgrl%r*ampl%r
                                end if
                                ampl = ampl * mp1_in%D0_norm

                            ! Stochastically coarse-grain the MP1 wavefunction by randomly rounding small amplitudes

                            ! Note attempt_to_spawn creates particles of opposite sign, hence set parent_sign = -1 
                            ! to undo this unwanted sign change.
                            ! pgen is always 1 as this is a deterministic calculation
                            ! TODO - cope with complex valued spawns
                            nspawn = attempt_to_spawn(rng, 1.0_p, spawn%cutoff, tijab%pop_real_factor, ampl%r, 1.0_p, -1_int_p)
                            ! cdet is not set, but since we're passing in f (the 'fexcit' argument), cdet is not required
                            if (nspawn /= 0_int_p) call create_spawned_particle(sys%basis, ref, cdet, excit, nspawn, 1, spawn, f)
                        end do
                    end do
                end do
            end do
            write (6,'(1X, A, ES17.10)') 'Deterministic MP2 correction energy: ', emp2
        end if


        call direct_annihilation(sys, rng, ref, annihilation_flags, tijab, spawn)
        if (parent) write (6,'()')

        call dealloc_spawn_t(spawn)

        ! Return sys in an unaltered state.
        call copy_sys_spin_info(sys_bak, sys)

    end subroutine sample_mp1_wfn

    subroutine mp1_in_t_json(js, mp1_in, terminal)

        ! Serialise a mp1_in_t object in JSON format.

        ! In/Out:
        !   js: json_out_t controlling the output unit and handling JSON internal state.  Unchanged on output.
        ! In:
        !   mp1_in: mp1_in_t object containing MP1 sampling input values (including any defaults set).
        !   terminal (optional): if true, this is the last entry in the enclosing JSON object.  Default: false.

        use json_out

        type(json_out_t), intent(inout) :: js
        type(mp1_in_t), intent(in) :: mp1_in
        logical, intent(in), optional :: terminal

        call json_object_init(js, 'mp1')
        call json_write_key(js, 'D0_norm', mp1_in%D0_norm)
        call json_write_key(js, 'state_size', mp1_in%state_size)
        call json_write_key(js, 'real_amplitudes', mp1_in%real_amplitudes)
        call json_write_key(js, 'spawn_cutoff', mp1_in%spawn_cutoff, terminal=.true.)
        call json_object_end(js, terminal)

    end subroutine mp1_in_t_json

end module mp1
