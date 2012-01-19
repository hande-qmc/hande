module fciqmc_common

! Module containing routines common to different fciqmc algorithms.

use fciqmc_data

implicit none

contains

    subroutine init_fciqmc()

        ! Initialisation for fciqmc calculations.
        ! Setup the spin polarisation for the system, initialise the RNG,
        ! allocate the required memory for the list of walkers and set the
        ! initial walker.

        use checking, only: check_allocate, check_deallocate
        use errors, only: stop_all
        use hashing, only: murmurhash_bit_string
        use parallel, only: iproc, nprocs, parent
        use proc_pointers
        use utils, only: int_fmt

        use annihilation, only: annihilate_main_list, annihilate_spawned_list, &
                                annihilate_main_list_initiator,                &
                                annihilate_spawned_list_initiator
        use basis, only: nbasis, basis_length, basis_fns, write_basis_fn, bit_lookup
        use basis, only: nbasis, basis_length, total_basis_length, basis_fns, write_basis_fn, basis_lookup, bit_lookup
        use calc, only: sym_in, ms_in, initiator_fciqmc, hfs_fciqmc_calc, ct_fciqmc_calc
        use calc, only: dmqmc_calc, doing_calc, doing_dmqmc_calc, dmqmc_energy, dmqmc_staggered_magnetisation
        use calc, only: dmqmc_energy_squared, dmqmc_correlation
        use dmqmc_procedures, only: init_dmqmc
        use determinants, only: encode_det, set_spin_polarisation, write_det
        use hamiltonian, only: get_hmatel
        use fciqmc_restart, only: read_restart
        use system, only: nel, nsites, ndim, system_type, hub_real, hub_k, heisenberg, staggered_magnetic_field
        use system, only: trial_function, neel_singlet, single_basis
        use symmetry, only: symmetry_orb_list
        use momentum_symmetry, only: gamma_sym, sym_table
        use utils, only: factorial_combination_1

        integer :: ierr
        integer :: i, D0_inv_proc, ipos, basis_find, occ_list0_inv(nel)
        integer :: step, size_main_walker, size_spawned_walker, nwalker_int, nwalker_real
        integer :: ref_sym ! the symmetry of the reference determinant
        integer(i0) :: f0_inv(basis_length)
        integer :: bit_position, bit_element

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
        if (doing_calc(initiator_fciqmc)) then
            spawned_size = spawned_size + 1
            spawned_parent = spawned_size
            annihilate_main_list_ptr => annihilate_main_list_initiator
            annihilate_spawned_list_ptr => annihilate_spawned_list_initiator
        else
            annihilate_main_list_ptr => annihilate_main_list
            annihilate_spawned_list_ptr => annihilate_spawned_list
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
            ! kinetic energy which satisfy the spin polarisation.
            ! Note: this is for testing only!  The symmetry input is currently
            ! ignored.
            call set_reference_det()

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
                write (6,'(a,64i4)') 'inv', occ_list0_inv
                    call encode_det(occ_list0_inv, f0_inv)
                    ! TODO: this needs to be checked for UHF upon merging
                    ! development branches.
!                    if (uhf) call stop_all('init_fciqmc','Check inversion for UHF systems.')
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
                    walker_population(:,tot_walkers) = 0.
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
            if (doing_calc(initiator_fciqmc)) then
                write (6,'(1X,a24)') 'Initiator method in use.'
                write (6,'(1X,a36,1X,"(",'//int_fmt(initiator_CAS(1),0)//',",",'//int_fmt(initiator_CAS(2),0)//'")")')  &
                    'CAS space of initiator determinants:',initiator_CAS
                write (6,'(1X,a66,'//int_fmt(initiator_population,1)//',/)') &
                    'Population for a determinant outside CAS space to be an initiator:', initiator_population
            end if
            write (6,'(1X,a49,/)') 'Information printed out every FCIQMC report loop:'
            write (6,'(1X,a66)') 'Instant shift: the shift calculated at the end of the report loop.'
            write (6,'(1X,a88)') 'Average shift: the running average &
                                                    &of the shift from when the shift was allowed to vary.'
            if (.not. doing_calc(dmqmc_calc)) then
                write (6,'(1X,a98)') 'Proj. Energy: projected energy averaged over the report loop. &
                                     &Calculated at the end of each cycle.'
                write (6,'(1X,a53)') 'Av. Proj. E: running average of the projected energy.'
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
                    write (6, '(1X,a106)') '\sum\rho_{ij}H_{ji}: The numerator of the estimator for the expectation &
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

    end subroutine init_fciqmc

    subroutine select_ref_det()

        ! Change the reference determinant to be the determinant with the
        ! greatest population if it exceeds some threshold relative to the
        ! current reference determinant.

        ! TODO: make compatible with HFS.

        use basis, only: basis_length
        use determinants, only: decode_det, write_det

        use parallel

        integer, parameter :: particle_type = 1
        integer :: in_data(2), out_data(2), i, fmax(basis_length), max_pop, ierr
        real(p) :: H00_max
        logical :: updated

        updated = .false.
        ! Find determinant with largest population.
        max_pop = 0
        do i = 1, tot_walkers
            if (abs(walker_population(particle_type,i)) > abs(max_pop)) then
                max_pop = walker_population(particle_type,i)
                fmax = walker_dets(:,i)
                H00_max = walker_data(particle_type, i)
            end if
        end do

        ! Only change reference determinant if the population is larger than the
        ! reference determinant by a given factor to avoid switching
        ! continuously between degenerate determinants.

        ! Note we don't broadcast the population of the new reference det as
        ! that is reset at the start of the next report loop anyway (and this
        ! routine should only be called at the end of the report loop).

#ifdef PARALLEL

        if (abs(max_pop) > ref_det_factor*abs(D0_population)) then
            in_data = (/ max_pop, iproc /)
        else if (iproc == D0_proc) then
            ! Ensure that D0_proc has the correct (average) population.
            in_data = (/ nint(D0_population), iproc /)
        else
            ! No det with sufficient population to become reference det on this
            ! processor.
            in_data = (/ 0, iproc /)
        end if

        call mpi_allreduce(in_data, out_data, 1, MPI_2INTEGER, MPI_MAXLOC, MPI_COMM_WORLD, ierr)

        if (out_data(1) /= nint(D0_population) .and. all(fmax /= f0)) then
            write (6,*) 'MPI_MAXLOC', out_data
            max_pop = out_data(1)
            updated = .true.
            D0_proc = out_data(2)
            f0 = fmax
            H00 = H00_max
            ! Broadcast updated data
            call mpi_bcast(f0, basis_length, mpi_det_integer, D0_proc, MPI_COMM_WORLD, ierr)
            call mpi_bcast(H00, 1, mpi_preal, D0_proc, MPI_COMM_WORLD, ierr)
        end if

#else

        if (abs(max_pop) > ref_det_factor*abs(D0_population) .and. all(fmax /= f0)) then
            updated = .true.
            f0 = fmax
            H00 = H00_max
        end if

#endif

        if (updated) then
            ! Update occ_list.
            call decode_det(f0, occ_list0)
            if (parent) then
                write (6,'(1X,"#",1X,62("-"))')
                write (6,'(1X,"#",1X,"Changed reference det to:",1X)',advance='no')
                call write_det(f0, new_line=.true.)
                write (6,'(1X,"#",1X,"Population on old reference det (averaged over report loop):",f10.2)') D0_population
                write (6,'(1X,"#",1X,"Population on new reference det:",27X,i8)') max_pop
                write (6,'(1X,"#",1X,"Care should be taken with accumulating statistics before this point.")')
                write (6,'(1X,"#",1X,62("-"))')
            end if
        end if

    end subroutine select_ref_det

    subroutine find_single_double_prob(occ_list, psingle, pdouble)

        ! Calculate the probabilities of selecting a single or double excitation
        ! from a given determinant.  We assume all possible excitations (i.e.
        ! those with Hamiltonian matrix elements which are not zero by symmetry)
        ! are equally likely, so this amounts to finding the number of possible
        ! (symmetry-allowed) single and double excitations.
        !
        ! In:
        !    occ_list: integer list of occupied spin-orbitals in a determinant, D.
        ! Out:
        !    psingle: probability of attempting to spawn on a determinant
        !             connected to D by a single excitation.
        !    pdouble: probability of attempting to spawn on a determinant
        !             connected to D by a double excitation.

        use basis, only: basis_fns
        use system, only: nel, system_type, hub_k, hub_real, heisenberg, read_in, sym0, nsym
        use point_group_symmetry, only: cross_product_pg_basis, cross_product_pg_sym, nbasis_sym_spin

        integer, intent(in) :: occ_list(nel)
        real(p), intent(out) :: psingle, pdouble

        integer :: i, j, virt_syms(2, sym0:nsym+(sym0-1)), nsingles, ndoubles, isyma, isymb, ims1, ims2

        select case(system_type)
        case(hub_k)
            ! Only double excitations
            psingle = 0.0_p
            pdouble = 1.0_p
        case(hub_real,heisenberg)
            ! Only single excitations
            psingle = 1.0_p
            pdouble = 0.0_p
        case(read_in)

            ! Count number of basis functions in each symmetry.
            virt_syms = nbasis_sym_spin
            do i = 1, nel
                ! Convert -1->1 and 1->2 for spin index in arrays.
                ims1 = (basis_fns(occ_list(i))%ms+3)/2
                virt_syms(ims1,basis_fns(occ_list(i))%sym) = virt_syms(ims1,basis_fns(occ_list(i))%sym) - 1
            end do

            ! Count number of possible single excitations from the supplied
            ! determinant.
            ! Symmetry and spin must be conserved.
            nsingles = 0
            do i = 1, nel
                ! Convert -1->1 and 1->2 for spin index in arrays.
                ims1 = (basis_fns(occ_list(i))%ms+3)/2
                ! Can't excite into already occupied orbitals.
                nsingles = nsingles + virt_syms(ims1,basis_fns(occ_list(i))%sym)
            end do

            ! Count number of possible double excitations from the supplied
            ! determinant.
            ndoubles = 0
            do i = 1, nel
                ! Convert -1->1 and 1->2 for spin index in arrays.
                ims1 = (basis_fns(occ_list(i))%ms+3)/2
                do j = i+1, nel
                    ! Convert -1->1 and 1->2 for spin index in arrays.
                    ims2 = (basis_fns(occ_list(j))%ms+3)/2
                    do isyma = sym0, nsym+(sym0-1)
                        ! Symmetry of the final orbital is determined (for Abelian
                        ! symmetries) from the symmetry of the first three.
                        isymb = cross_product_pg_sym(isyma, cross_product_pg_basis(occ_list(i),occ_list(j)))
                        if (isyma == isymb) then
                            if (ims1 == ims2) then
                                ! Cannot excit 2 electrons into the same spin-orbital.
                                ! Need to avoid double counting.
                                !  => number of unique pairs is identical to
                                !  number of elements in the strictly lower
                                !  triangle of a square matrix.
                                ndoubles = ndoubles + (virt_syms(ims1,isyma)*(virt_syms(ims2,isymb)-1))/2
                            else
                                ndoubles = ndoubles + virt_syms(ims1,isyma)*virt_syms(ims2,isymb)
                            end if
                        else if (isyma < isymb) then
                            ! isyma < isymb to avoid double counting.
                            ndoubles = ndoubles + virt_syms(ims1,isyma)*virt_syms(ims2,isymb)
                            ! can also have the opposite spin structure of
                            ! occupied orbitals have different spins.
                            if (ims1 /= ims2) ndoubles = ndoubles + virt_syms(ims2,isyma)*virt_syms(ims1,isymb)
                        end if
                    end do
                end do
            end do

            psingle = real(nsingles,p)/(nsingles+ndoubles)
            pdouble = real(ndoubles,p)/(nsingles+ndoubles)

        end select

    end subroutine find_single_double_prob

    subroutine initial_fciqmc_status()

        ! Calculate the projected energy based upon the initial walker
        ! distribution (either via a restart or as set during initialisation)
        ! and print out.

        use parallel
        use proc_pointers, only: update_proj_energy_ptr
        integer :: idet
        integer :: ntot_particles
#ifdef PARALLEL
        integer :: ierr
        real(p) :: proj_energy_sum
#endif

        ! Calculate the projected energy based upon the initial walker
        ! distribution.  proj_energy and D0_population are both accumulated in
        ! update_proj_energy.
        proj_energy = 0.0_p
        D0_population = 0
        do idet = 1, tot_walkers
            call update_proj_energy_ptr(idet)
        end do

#ifdef PARALLEL
        call mpi_allreduce(proj_energy, proj_energy_sum, 1, mpi_preal, MPI_SUM, MPI_COMM_WORLD, ierr)
        proj_energy = proj_energy_sum
        call mpi_allreduce(nparticles, ntot_particles, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
        ntot_particles = nparticles(1)
#endif

        proj_energy = proj_energy/D0_population

        if (parent) then
            ! See also the format used in write_fciqmc_report if this is changed.
            ! We prepend a # to make it easy to skip this point when do data
            ! analysis.
            write (6,'(1X,"#",3X,i8,2X,4(es17.10,2X),es17.10,4X,i11,6X,a3,3X,a3)') &
                    mc_cycles_done, shift, 0.0_p, proj_energy, 0.0_p, D0_population, ntot_particles,'n/a','n/a'
        end if

    end subroutine initial_fciqmc_status

    subroutine load_balancing_report()

        ! Print out a load-balancing report when run in parallel showing how
        ! determinants and walkers/particles are distributed over the processors.

#ifdef PARALLEL
        use annihilation, only: annihilation_comms_time
        use parallel

        integer(lint) :: load_data(nprocs)
        integer :: ierr
        real(dp) :: comms(nprocs)

        if (nprocs > 1) then
            if (parent) then
                write (6,'(1X,a14,/,1X,14("^"),/)') 'Load balancing'
                write (6,'(1X,a77,/)') "The final distribution of walkers and determinants across the processors was:"
            endif
            call mpi_gather(nparticles, 1, mpi_integer8, load_data, 1, mpi_integer8, 0, MPI_COMM_WORLD, ierr)
            if (parent) then
                write (6,'(1X,a34,6X,i8)') 'Min # of particles on a processor:', minval(load_data)
                write (6,'(1X,a34,6X,i8)') 'Max # of particles on a processor:', maxval(load_data)
                write (6,'(1X,a35,5X,f11.2)') 'Mean # of particles on a processor:', real(sum(load_data), p)/nprocs
            end if
            call mpi_gather(tot_walkers, 1, mpi_integer, load_data, 1, mpi_integer8, 0, MPI_COMM_WORLD, ierr)
            call mpi_gather(annihilation_comms_time, 1, mpi_real8, comms, 1, mpi_real8, 0, MPI_COMM_WORLD, ierr)
            if (parent) then
                write (6,'(1X,a37,3X,i8)') 'Min # of determinants on a processor:', minval(load_data)
                write (6,'(1X,a37,3X,i8)') 'Max # of determinants on a processor:', maxval(load_data)
                write (6,'(1X,a38,2X,f11.2)') 'Mean # of determinants on a processor:', real(sum(load_data), p)/nprocs
                write (6,'()')
                write (6,'(1X,a38,5X,f8.2,a1)') 'Min time take by walker communication:', minval(comms),'s'
                write (6,'(1X,a38,5X,f8.2,a1)') 'Max time take by walker communication:', maxval(comms),'s'
                write (6,'(1X,a39,4X,f8.2,a1)') 'Mean time take by walker communication:', real(sum(comms), p)/nprocs,'s'
                write (6,'()')
            end if
        end if
#endif

    end subroutine load_balancing_report

end module fciqmc_common
