module qmc

! Launcher and initialisation routines for the various QMC algorithms.

use fciqmc_data

implicit none

contains

! --- QMC wrapper ---

    subroutine do_qmc()

        use calc

        use ccmc, only: do_ccmc
        use ct_fciqmc, only: do_ct_fciqmc
        use dmqmc, only: do_dmqmc
        use fciqmc, only: do_fciqmc
        use folded_spectrum_utils, only: init_folded_spectrum
        use ifciqmc, only: init_ifciqmc
        use hellmann_feynman_sampling, only: init_hellmann_feynman_sampling, do_hfs_fciqmc
        use energy_evaluation, only: update_proj_hfs_hub_k

        real(dp) :: hub_matel

        ! Initialise data
        call init_qmc()

        ! Initialise procedure pointers
        call init_proc_pointers()

        ! Calculation-specifc initialisation and then run QMC calculation.

        if (initiator_approximation) call init_ifciqmc()

        if (doing_calc(dmqmc_calc)) then
            call do_dmqmc()
        else if (doing_calc(ct_fciqmc_calc)) then
            call do_ct_fciqmc(hub_matel)
        else if (doing_calc(ccmc_calc)) then
            call do_ccmc()
        else
            ! Doing FCIQMC calculation (of some sort) using the original
            ! timestep algorithm.
            if (doing_calc(folded_spectrum)) call init_folded_spectrum()
            if (doing_calc(hfs_fciqmc_calc)) then
                call init_hellmann_feynman_sampling()
                call do_hfs_fciqmc(update_proj_hfs_hub_k)
            else
                call do_fciqmc()
            end if
        end if

    end subroutine do_qmc

! --- Initialisation routines ---

    subroutine init_qmc()

        ! Initialisation for fciqmc calculations.
        ! Setup the spin polarisation for the system, initialise the RNG,
        ! allocate the required memory for the list of walkers and set the
        ! initial walker.

        use checking, only: check_allocate, check_deallocate
        use errors, only: stop_all
        use hashing, only: murmurhash_bit_string
        use parallel, only: iproc, nprocs, parent
        use utils, only: int_fmt

        use annihilation, only: annihilate_main_list, annihilate_spawned_list, &
                                annihilate_main_list_initiator,                &
                                annihilate_spawned_list_initiator
        use basis, only: nbasis, basis_length, basis_fns, write_basis_fn, bit_lookup
        use basis, only: nbasis, basis_length, total_basis_length, basis_fns, write_basis_fn, basis_lookup, bit_lookup
        use calc, only: sym_in, ms_in, initiator_approximation, fciqmc_calc, hfs_fciqmc_calc, ct_fciqmc_calc
        use calc, only: dmqmc_calc, doing_calc, doing_dmqmc_calc, dmqmc_energy, dmqmc_staggered_magnetisation
        use calc, only: dmqmc_energy_squared, dmqmc_correlation
        use dmqmc_procedures, only: init_dmqmc
        use determinants, only: encode_det, set_spin_polarisation, write_det
        use hamiltonian, only: get_hmatel
        use qmc_common, only: find_single_double_prob
        use fciqmc_restart, only: read_restart
        use reference_determinant, only: set_reference_det
        use system, only: nel, nsites, ndim, system_type, hub_real, hub_k, heisenberg, staggered_magnetic_field
        use system, only: trial_function, neel_singlet, single_basis, sym_max
        use symmetry, only: symmetry_orb_list
        use momentum_symmetry, only: gamma_sym, sym_table
        use utils, only: factorial_combination_1

        integer :: ierr
        integer :: i, D0_inv_proc, ipos, occ_list0_inv(nel)
        integer :: step, size_main_walker, size_spawned_walker, nwalker_int, nwalker_real
        integer :: ref_sym ! the symmetry of the reference determinant
        integer(i0) :: f0_inv(basis_length)

        if (parent) write (6,'(1X,a6,/,1X,6("-"),/)') 'FCIQMC'

        ! Array sizes depending upon FCIQMC algorithm.
        sampling_size = 1
        spawned_size = total_basis_length + 1
        spawned_pop = spawned_size
        if (doing_calc(hfs_fciqmc_calc)) then
            spawned_size = spawned_size + 1
            spawned_hf_pop = spawned_size
            sampling_size = sampling_size + 1
        else
            spawned_hf_pop = spawned_size
        end if
        if (initiator_approximation) then
            spawned_size = spawned_size + 1
            spawned_parent = spawned_size
        end if

        ! Each determinant occupies basis_length kind=i0 integers,
        ! sampling_size integers and sampling_size kind=p reals.
        ! If the Neel singlet state is used as the reference state for the
        ! projected estimator, then a further 2 reals are used per
        ! determinant.
        if (trial_function == neel_singlet) then
            info_size = 2
        else
            info_size = 0
        end if
        nwalker_int = sampling_size
        nwalker_real = sampling_size + info_size

        ! Thus the number of bits occupied by each determinant in the main
        ! walker list is given by basis_length*i0_length+nwalker_int*32+nwalker_real*32
        ! (*64 if double precision).  The number of bytes is simply 1/8 this.
#ifdef SINGLE_PRECISION
        size_main_walker = total_basis_length*i0_length/8 + nwalker_int*4 + nwalker_real*4
#else
        size_main_walker = total_basis_length*i0_length/8 + nwalker_int*4 + nwalker_real*8
#endif
        if (walker_length < 0) then
            ! Given in MB.  Convert.  Note: important to avoid overflow in the
            ! conversion!
            walker_length = int((-real(walker_length,p)*10**6)/size_main_walker)
        end if

        ! Each spawned_walker occupies spawned_size kind=i0 integers.
        size_spawned_walker = spawned_size*i0_length/8
        if (spawned_walker_length < 0) then
            ! Given in MB.  Convert.
            ! Note that we store 2 arrays.
            spawned_walker_length = int((-real(spawned_walker_length,p)*10**6)/(2*size_spawned_walker))
        end if

        if (parent) then
            write (6,'(1X,a53,f7.2)') &
                'Memory allocated per core for main walker list (MB): ', &
                size_main_walker*real(walker_length,p)/10**6
            write (6,'(1X,a57,f7.2)') &
                'Memory allocated per core for spawned walker lists (MB): ', &
                size_spawned_walker*real(2*spawned_walker_length,p)/10**6
            write (6,'(1X,a48,'//int_fmt(walker_length,1)//')') &
                'Number of elements per core in main walker list:', walker_length
            write (6,'(1X,a51,'//int_fmt(spawned_walker_length,1)//',/)') &
                'Number of elements per core in spawned walker list:', spawned_walker_length
        end if

        ! Allocate main walker lists.
        allocate(nparticles(sampling_size), stat=ierr)
        call check_allocate('nparticles', sampling_size, ierr)
        allocate(walker_dets(total_basis_length,walker_length), stat=ierr)
        call check_allocate('walker_dets', basis_length*walker_length, ierr)
        allocate(walker_population(sampling_size,walker_length), stat=ierr)
        call check_allocate('walker_population', sampling_size*walker_length, ierr)
        allocate(walker_data(sampling_size+info_size,walker_length), stat=ierr)
        call check_allocate('walker_data', size(walker_data), ierr)

        ! Allocate spawned walker lists and spawned walker times (ct_fciqmc only)
        if (mod(spawned_walker_length, nprocs) /= 0) then
            if (parent) write (6,'(1X,a68)') 'spawned_walker_length is not a multiple of the number of processors.'
            spawned_walker_length = ceiling(real(spawned_walker_length)/nprocs)*nprocs
            if (parent) write (6,'(1X,a35,'//int_fmt(spawned_walker_length,1)//',a1,/)') &
                                        'Increasing spawned_walker_length to',spawned_walker_length,'.'
        end if
        allocate(spawned_walkers1(spawned_size,spawned_walker_length), stat=ierr)
        call check_allocate('spawned_walkers1',spawned_size*spawned_walker_length,ierr)
        spawned_walkers => spawned_walkers1
        if (doing_calc(ct_fciqmc_calc)) then
            allocate(spawn_times(spawned_walker_length),stat=ierr)
            call check_allocate('spawn_times',spawned_walker_length,ierr)
        end if
        ! Allocate scratch space for doing communication.
        allocate(spawned_walkers2(spawned_size,spawned_walker_length), stat=ierr)
        call check_allocate('spawned_walkers2',spawned_size*spawned_walker_length,ierr)
        spawned_walkers_recvd => spawned_walkers2

        ! Set spawning_head to be the same size as spawning_block_start.
        allocate(spawning_head(0:max(1,nprocs-1)), stat=ierr)
        call check_allocate('spawning_head',max(2,nprocs),ierr)

        ! Find the start position within the spawned walker lists for each
        ! processor.
        ! spawning_block_start(1) should contain the number of elements allocated
        ! for each processor so we allow it to be accessible even if the number
        ! of processors is 1.
        allocate(spawning_block_start(0:max(1,nprocs-1)), stat=ierr)
        call check_allocate('spawning_block_start',max(2,nprocs),ierr)
        step = spawned_walker_length/nprocs
        forall (i=0:nprocs-1) spawning_block_start(i) = i*step

        ! Set spin variables for non-Heisenberg systems
        if (system_type /= heisenberg) call set_spin_polarisation(ms_in)

        ! Set initial walker population.
        ! occ_list could be set and allocated in the input.
        allocate(f0(basis_length), stat=ierr)
        call check_allocate('f0',basis_length,ierr)
        if (restart) then
            if (.not.allocated(occ_list0)) then
                allocate(occ_list0(nel), stat=ierr)
                call check_allocate('occ_list0',nel,ierr)
            end if
            call read_restart()
        else

            ! Reference det

            ! Set the reference determinant to be the spin-orbitals with the lowest
            ! single-particle eigenvalues which satisfy the spin polarisation.
            ! Note: this is for testing only!  The symmetry input is currently
            ! ignored.
            if (sym_in < sym_max) then
                call set_reference_det(occ_list0, .false., sym_in)
            else
                call set_reference_det(occ_list0, .false.)
            end if

            call encode_det(occ_list0, f0)

            ! In general FCIQMC, we start with psips only on the
            ! reference determinant, so set tot_walkers = 1 and
            ! initialise walker_population. For DMQMC, this is
            ! not required, as psips are spawned along the diagonal
            ! initially.
            if (.not.doing_calc(dmqmc_calc)) then
                tot_walkers = 1
                ! Zero all populations...
                walker_population(:,tot_walkers) = 0
                ! Set initial population of Hamiltonian walkers.
                walker_population(1,tot_walkers) = nint(D0_population)
                ! Set the bitstring of this psip to be that of the
                ! reference state.
                walker_dets(:,tot_walkers) = f0
            end if

            ! Energy of reference determinant.
            H00 = get_hmatel(f0,f0)

            ! Determine and set properties for the reference state which we start on.
            ! (For DMQMC, we do not start on the reference state, and so this is not
            ! required. Psips are initialised along the diagonal in DMQMC. See
            ! dmqmc_procedures).
            if (.not.doing_calc(dmqmc_calc)) then
                ! By definition, when using a single determinant as a reference state:
                walker_data(1,tot_walkers) = 0.0_p
                ! Or if not using a single determinant:
                if (trial_function == neel_singlet) then
                    ! Set the Neel state data for the reference state, if it is being used.
                    walker_data(1,tot_walkers) = H00
                    H00 = 0.0_p
                    walker_data(sampling_size+1,tot_walkers) = nsites/2
                    ! For a rectangular bipartite lattice, nbonds = ndim*nsites.
                    ! The Neel state cannot be used for non-bipartite lattices.
                    walker_data(sampling_size+2,tot_walkers) = ndim*nsites
                end if

                ! Finally, we need to check if the reference determinant actually
                ! belongs on this processor.
                ! If it doesn't, set the walkers array to be empty.
                if (nprocs > 1) then
                    D0_proc = modulo(murmurhash_bit_string(f0, basis_length), nprocs)
                    if (D0_proc /= iproc) tot_walkers = 0
                else
                    D0_proc = iproc
                end if
            end if

            ! For the Heisenberg model and open shell systems, it is often useful to
            ! have psips start on both the reference state and the spin-flipped version.
            if (init_spin_inv_D0) then

                ! Need to handle the Heisenberg model (consisting of spinors on
                ! lattice sites) and electron systems differently, as the
                ! Heisenberg model has no concept of unoccupied basis
                ! functions/holes.
                select case (system_type)
                case (heisenberg)
                    ! Flip all spins in f0 to get f0_inv
                    f0_inv = not(f0)
                    ! In general, the basis bit string has some padding at the
                    ! end which must be unset.  We need to clear this...
                    ! Loop over all bits after the last basis function.
                    i = bit_lookup(2,nbasis)
                    do ipos = bit_lookup(1,nbasis)+1, i0_end
                        f0_inv(i) = ibclr(f0_inv(i), ipos)
                    end do
                case default
                    ! Swap each basis function for its spin-inverse
                    ! This looks somewhat odd, but relies upon basis
                    ! functions alternating down (Ms=-1) and up (Ms=1).
                    do i = 1, nel
                        if (mod(occ_list0(i),2) == 1) then
                            occ_list0_inv(i) = occ_list0(i) + 1
                        else
                            occ_list0_inv(i) = occ_list0(i) - 1
                        end if
                    end do
                    call encode_det(occ_list0_inv, f0_inv)
                end select

                if (nprocs > 1) then
                    D0_inv_proc = modulo(murmurhash_bit_string(f0_inv, basis_length), nprocs)
                else
                    D0_inv_proc = iproc
                end if

                ! Store if not identical to reference det.
                if (D0_inv_proc == iproc .and. any(f0 /= f0_inv)) then
                    tot_walkers = tot_walkers + 1
                    ! Zero all populations for this determinant.
                    walker_population(:,tot_walkers) = 0
                    ! Set the population for this basis function.
                    walker_population(1,tot_walkers) = nint(D0_population)
                    walker_data(1,tot_walkers) = get_hmatel(f0,f0) - H00
                    select case(system_type)
                    case(heisenberg)
                        if (trial_function /= single_basis) then
                            walker_data(1,tot_walkers) = 0
                        else
                            walker_data(1,tot_walkers) = get_hmatel(f0,f0) - H00
                        end if
                    case default
                        walker_data(1,tot_walkers) = get_hmatel(f0,f0) - H00
                    end select
                    walker_dets(:,tot_walkers) = f0_inv
                    ! If we are using the Neel state as a reference in the
                    ! Heisenberg model, then set the required data.
                    if (trial_function == neel_singlet) then
                        walker_data(sampling_size+1,tot_walkers) = 0
                        walker_data(sampling_size+2,tot_walkers) = 0
                    end if
                end if
            end if

        end if ! End of initialisation of reference state(s)/restarting from previous calculations

        ! Total number of particles on processor.
        ! Probably should be handled more simply by setting it to be either 0 or
        ! D0_population or obtaining it from the restart file, as appropriate.
        forall (i=1:sampling_size) nparticles(i) = sum(abs(walker_population(i,:tot_walkers)))

        ! calculate the reference determinant symmetry
        ref_sym = symmetry_orb_list(occ_list0)

        ! If not set at input, set probability of selecting single or double
        ! excitations based upon the reference determinant and assume other
        ! determinants have a roughly similar ratio of single:double
        ! excitations.
        if (pattempt_single < 0 .or. pattempt_double < 0) then
            call find_single_double_prob(occ_list0, pattempt_single, pattempt_double)
        else
            ! renormalise just in case input wasn't
            pattempt_single = pattempt_single/(pattempt_single+pattempt_double)
            pattempt_double = 1.0_p - pattempt_single
        end if

        ! Calculate all the possible different amplitudes for the Neel singlet state
        ! and store them in an array
        if (trial_function == neel_singlet) then
            allocate(neel_singlet_amp(-1:(nsites/2)+1), stat=ierr)
            call check_allocate('neel_singlet_amp',(nsites/2)+1,ierr)

            neel_singlet_amp(-1) = 0
            neel_singlet_amp((nsites/2)+1) = 0
            do i=0,(nsites/2)
                neel_singlet_amp(i) = factorial_combination_1( (nsites/2)-i , i )
                neel_singlet_amp(i) = -(2*mod(i,2)-1) * neel_singlet_amp(i)
            end do
        end if

        ! When doing a DMQMC calculation, call a routine to initialise all the required
        ! arrays, ie to store thermal quantities, and to initalise reduced density matrix
        ! quantities if necessary.
        if (doing_calc(dmqmc_calc)) then
            call init_dmqmc()
        end if

        if (parent) then
            write (6,'(1X,a29,1X)',advance='no') 'Reference determinant, |D0> ='
            call write_det(f0, new_line=.true.)
            write (6,'(1X,a16,f20.12)') 'E0 = <D0|H|D0> =',H00
            write(6,'(1X,a34)',advance='no') 'Symmetry of reference determinant:'
            select case(system_type)
            case (hub_k)
                call write_basis_fn(basis_fns(2*ref_sym), new_line=.true., print_full=.false.)
            case default
                write(6,'(i2)') ref_sym
            end select
            write (6,'(1X,a46,1X,f8.4)') 'Probability of attempting a single excitation:', pattempt_single
            write (6,'(1X,a46,1X,f8.4)') 'Probability of attempting a double excitation:', pattempt_double
            write (6,'(1X,a44,1X,f11.4,/)') &
                              'Initial population on reference determinant:',D0_population
            write (6,'(1X,a68,/)') 'Note that FCIQMC calculates the correlation energy relative to |D0>.'
            if (initiator_approximation) then
                write (6,'(1X,a24)') 'Initiator method in use.'
                if (doing_calc(fciqmc_calc)) &
                    write (6,'(1X,a36,1X,"(",'//int_fmt(initiator_CAS(1),0)//',",",'//int_fmt(initiator_CAS(2),0)//'")")')  &
                    'CAS space of initiator determinants:',initiator_CAS
                write (6,'(1X,a66,'//int_fmt(initiator_population,1)//',/)') &
                    'Population for a determinant outside CAS space to be an initiator:', initiator_population
            end if
            write (6,'(1X,a49,/)') 'Information printed out every FCIQMC report loop:'
            write (6,'(1X,a66)') 'Instant shift: the shift calculated at the end of the report loop.'
            if (.not. doing_calc(dmqmc_calc)) then
                write (6,'(1X,a98)') 'Proj. Energy: projected energy averaged over the report loop. &
                                     &Calculated at the end of each cycle.'
                write (6,'(1X,a54)') '# D0: current population at the reference determinant.'
            else
                write (6, '(1X,a83)') 'Trace: The current total population on the diagonal elements of the &
                                     &density matrix.'
                if (doing_dmqmc_calc(dmqmc_energy)) then
                    write (6, '(1X,a92)') '\sum\rho_{ij}H_{ji}: The numerator of the estimator for the expectation &
                                         &value of the energy.'
                end if
                if (doing_dmqmc_calc(dmqmc_energy_squared)) then
                    write (6, '(1X,a100)') '\sum\rho_{ij}H2{ji}: The numerator of the estimator for the expectation &
                                         &value of the energy squared.'
                end if
                if (doing_dmqmc_calc(dmqmc_correlation)) then
                    write (6, '(1X,a111)') '\sum\rho_{ij}H_{ji}: The numerator of the estimator for the expectation &
                                         &value of the spin correlation function.'
                end if
                if (doing_dmqmc_calc(dmqmc_staggered_magnetisation)) then
                    write (6, '(1X,a109)') '\sum\rho_{ij}M2{ji}: The numerator of the estimator for the expectation &
                                         &value of the staggered magnetisation.'
                end if
            end if
            write (6,'(1X,a49)') '# particles: current total population of walkers.'
            write (6,'(1X,a56,/)') 'R_spawn: average rate of spawning across all processors.'
        end if

    end subroutine init_qmc

    subroutine init_proc_pointers()

        ! Set function pointers for QMC calculations.

        ! System and calculation data
        use calc
        use system

        ! Procedures to be pointed to.
        use annihilation
        use death, only: stochastic_death
        use determinants
        use dmqmc_estimators
        use dmqmc_procedures
        use energy_evaluation
        use excit_gen_mol
        use excit_gen_hub_k
        use excit_gen_real_lattice
        use folded_spectrum_utils, only: fs_spawner, fs_stochastic_death
        use hamiltonian_chung_landau, only: slater_condon0_chung_landau
        use hamiltonian_hub_k, only: slater_condon0_hub_k
        use hamiltonian_hub_real, only: slater_condon0_hub_real
        use hamiltonian_heisenberg, only: diagonal_element_heisenberg, diagonal_element_heisenberg_staggered
        use hamiltonian_molecular, only: slater_condon0_mol
        use heisenberg_estimators
        use ifciqmc, only: set_parent_flag, set_parent_flag_dummy
        use importance_sampling
        use spawning

        ! Procedure pointers
        use proc_pointers

        ! Utilities
        use errors, only: stop_all

        ! 0. In general, use the default spawning routine.
        spawner_ptr => spawn

        ! 1. Set system-specific procedure pointers.
        !     * projected energy estimator
        !     * diagonal hamiltonian matrix element evaluation
        !     * spawning
        !     * excitation generators
        select case(system_type)
        case(hub_k)

            decoder_ptr => decode_det_spinocc_spinunocc
            update_proj_energy_ptr => update_proj_energy_hub_k
            sc0_ptr => slater_condon0_hub_k

            spawner_ptr => spawn_lattice_split_gen
            if (no_renorm) then
                gen_excit_ptr => gen_excit_hub_k_no_renorm
                gen_excit_init_ptr => gen_excit_init_hub_k_no_renorm
                gen_excit_finalise_ptr => gen_excit_finalise_hub_k_no_renorm
            else
                gen_excit_ptr => gen_excit_hub_k
                gen_excit_init_ptr => gen_excit_init_hub_k
                gen_excit_finalise_ptr => gen_excit_finalise_hub_k
            end if

        case(hub_real, chung_landau)

            ! The Hubbard model in a local orbital basis and the Chung--Landau
            ! Hamiltonian have the same off-diagonal operator so we can use the
            ! same excitation generators and just a different function for the
            ! diagonal matrix elements.
            ! Note that the Chung--Landau model (as in Phys Rev B 85 (2012)
            ! 115115) is contains spinless fermions.
            decoder_ptr => decode_det_occ
            update_proj_energy_ptr => update_proj_energy_hub_real
            if (system_type == hub_real) then
                sc0_ptr => slater_condon0_hub_real
            else
                sc0_ptr => slater_condon0_chung_landau
            end if

            if (no_renorm) then
                gen_excit_ptr => gen_excit_hub_real_no_renorm
            else
                gen_excit_ptr => gen_excit_hub_real
            end if

        case(heisenberg)

            ! Only need occupied orbitals list, as for the real Hubbard case.
            decoder_ptr => decode_det_occ
            ! Set which trial wavefunction to use for the energy estimator.
            select case(trial_function)
            case (single_basis)
                update_proj_energy_ptr => update_proj_energy_heisenberg_basic
            case (neel_singlet)
                update_proj_energy_ptr => update_proj_energy_heisenberg_neel_singlet
            end select

            ! Set whether the applied staggered magnetisation is non-zero.
            if (abs(staggered_magnetic_field) > 0.0_p) then
                sc0_ptr => diagonal_element_heisenberg_staggered
            else
                sc0_ptr => diagonal_element_heisenberg
            end if

            ! Set which guiding wavefunction to use, if requested.
            if (no_renorm) then
                gen_excit_ptr => gen_excit_heisenberg_no_renorm
            else
                    gen_excit_ptr => gen_excit_heisenberg
            end if
            select case(guiding_function)
            case (neel_singlet_guiding)
                spawner_ptr => spawn_importance_sampling
                trial_fn_ptr => neel_trial_state
            end select

        case(read_in)

            update_proj_energy_ptr => update_proj_energy_mol
            sc0_ptr => slater_condon0_mol

            if (no_renorm) then
                gen_excit_ptr => gen_excit_mol_no_renorm
                decoder_ptr => decode_det_occ
            else
                gen_excit_ptr => gen_excit_mol
                decoder_ptr => decode_det_occ_symunocc
            end if

        case(ueg)
            call stop_all('init_proc_pointers','QMC not implemented for the UEG yet.')
        case default
            call stop_all('init_proc_pointers','QMC not implemented for this system yet.')
        end select

        ! 2. Set calculation-specific procedure pointers

        ! a) initiator-approximation
        if (initiator_approximation) then
            set_parent_flag_ptr => set_parent_flag
            if (truncate_space) then
                create_spawned_particle_ptr => create_spawned_particle_initiator_truncated
            else
                create_spawned_particle_ptr => create_spawned_particle_initiator
            end if
            annihilate_main_list_ptr => annihilate_main_list_initiator
            annihilate_spawned_list_ptr => annihilate_spawned_list_initiator
        else
            set_parent_flag_ptr => set_parent_flag_dummy
            if (truncate_space) then
                create_spawned_particle_ptr => create_spawned_particle_truncated
            else
                create_spawned_particle_ptr => create_spawned_particle
            end if
            annihilate_main_list_ptr => annihilate_main_list
            annihilate_spawned_list_ptr => annihilate_spawned_list
        end if

        ! b) folded-spectrum
        if (doing_calc(folded_spectrum)) then
            spawner_ptr => fs_spawner
            death_ptr => fs_stochastic_death
        else
            death_ptr => stochastic_death
        end if

        ! c) density-matrix
        if (doing_calc(dmqmc_calc)) then

            ! Spawned particle creation. 
            if (truncate_space) then
                create_spawned_particle_dm_ptr => create_spawned_particle_truncated_density_matrix
            else
                create_spawned_particle_dm_ptr => create_spawned_particle_density_matrix
            end if

            ! Expectation values.
            select case(system_type)
            case(heisenberg)
                if (doing_dmqmc_calc(dmqmc_energy)) update_dmqmc_energy_ptr => dmqmc_energy_heisenberg
                if (doing_dmqmc_calc(dmqmc_energy_squared)) &
                                         update_dmqmc_energy_squared_ptr => dmqmc_energy_squared_heisenberg
                if (doing_dmqmc_calc(dmqmc_correlation)) update_dmqmc_correlation_ptr => &
                                     dmqmc_correlation_function_heisenberg
                if (doing_dmqmc_calc(dmqmc_staggered_magnetisation)) &
                                         update_dmqmc_stag_mag_ptr => dmqmc_stag_mag_heisenberg
            end select

        end if

    end subroutine init_proc_pointers

end module qmc
