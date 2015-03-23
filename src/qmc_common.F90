module qmc_common

! Module containing routines common to different qmc algorithms.

use fciqmc_data

implicit none

contains

! --- Utility routines ---

    subroutine select_ref_det(sys, ref_det_factor, reference, psip_list, qs)

        ! Change the reference determinant to be the determinant with the
        ! greatest population if it exceeds some threshold relative to the
        ! current reference determinant.

        ! Note this currently only looks at the Hamiltonian population.  The
        ! setting of multiple reference determinants (e.g. different references
        ! for Hamiltonian walkers and Hellmann-Feynmann walkers) is currently not
        ! supported.  It is not clear if there is such a need as a good
        ! reference determinant for the Hamiltonian space should also be
        ! important (if not crucial!) in the H-F space.

        ! In:
        !    sys: system being studied.
        !    ref_det_factor: factor by which population on a determinant must
        !        exceed the reference population in order to be accepted as
        !        the new reference.
        ! In/Out:
        !    reference: reference determinant.
        !    psip_list: particle_t object containing the main particle list.  On
        !       output the diagonal Hamiltonian matrix elements are updated with
        !       to be relative to the new reference determinant (if applicable).

        use calc, only: doing_calc, hfs_fciqmc_calc
        use determinants, only: decode_det, write_det
        use system, only: sys_t

        use parallel
        use errors, only: stop_all
        use qmc_data, only: reference_t, particle_t, qmc_state_t

        type(sys_t), intent(in) :: sys
        real(p), intent(in) :: ref_det_factor
        type(reference_t), intent(inout) :: reference
        type(particle_t), intent(inout) :: psip_list
        type(qmc_state_t), intent(inout) :: qs

        integer, parameter :: particle_type = 1
        integer :: i, D0_proc
        integer(i0), allocatable :: fmax(:)
        integer(int_p) :: max_pop
#ifdef PARALLEL
        integer(int_p) :: in_data(2), out_data(2)
        integer :: ierr
#endif
        real(p) :: H00_max, H00_old
        logical :: updated

        allocate(fmax(lbound(psip_list%states, dim=1):ubound(psip_list%states, dim=1)))

        H00_old = reference%H00

        updated = .false.
        ! Find determinant with largest population.
        max_pop = 0_int_p
        do i = 1, psip_list%nstates
            if (abs(psip_list%pops(particle_type,i)) > abs(max_pop)) then
                max_pop = psip_list%pops(particle_type,i)
                fmax = psip_list%states(:,i)
                H00_max = psip_list%dat(particle_type, i)
            end if
        end do

        ! Only change reference determinant if the population is larger than the
        ! reference determinant by a given factor to avoid switching
        ! continuously between degenerate determinants.

        ! Note we don't broadcast the population of the new reference det as
        ! that is reset at the start of the next report loop anyway (and this
        ! routine should only be called at the end of the report loop).

#ifdef PARALLEL

        if (all(fmax == reference%f0)) then
            ! Max population on this processor is already the reference.  Don't change.
            in_data = (/ 0_int_p, int(iproc,int_p) /)
        else if (abs(max_pop) > ref_det_factor*abs(qs%D0_population)) then
            in_data = (/ max_pop, int(iproc,int_p) /)
        else
            ! No det with sufficient population to become reference det on this
            ! processor.
            in_data = (/ 0_int_p, int(iproc, int_p) /)
        end if

        call mpi_allreduce(in_data, out_data, 2, mpi_pop_integer, MPI_MAXLOC, MPI_COMM_WORLD, ierr)

        if (out_data(1) /= 0) then
            max_pop = out_data(1)
            updated = .true.
            D0_proc = out_data(2)
            reference%f0 = fmax
            reference%H00 = H00_max
            ! Broadcast updated data
            call mpi_bcast(reference%f0, size(reference%f0), mpi_det_integer, D0_proc, MPI_COMM_WORLD, ierr)
            call mpi_bcast(reference%H00, 1, mpi_preal, D0_proc, MPI_COMM_WORLD, ierr)
        end if

#else

        if (abs(max_pop) > ref_det_factor*abs(qs%D0_population) .and. any(fmax /= reference%f0)) then
            updated = .true.
            reference%f0 = fmax
            reference%H00 = H00_max
        end if

#endif

        if (updated) then
            ! Update occ_list.
            call decode_det(sys%basis, reference%f0, reference%occ_list0)
            ! psip_list%dat(1,i) holds <D_i|H|D_i> - H00_old.  Update.
            ! H00 is currently <D_0|H|D_0> - H00_old.
            ! Want psip_list%dat(1,i) to be <D_i|H|D_i> - <D_0|H|D_0>
            ! We'll fix H00 later and avoid an extra nstates*additions.
            do i = 1, psip_list%nstates
                psip_list%dat(1,i) = psip_list%dat(1,i) - reference%H00
            end do
            ! Now set H00 = <D_0|H|D_0> so that future references to it are
            ! correct.
            reference%H00 = reference%H00 + H00_old
            if (doing_calc(hfs_fciqmc_calc)) call stop_all('select_ref_det', 'Not implemented for HFS.')
            if (parent) then
                write (6,'(1X,"#",1X,62("-"))')
                write (6,'(1X,"#",1X,"Changed reference det to:",1X)',advance='no')
                call write_det(sys%basis, sys%nel, reference%f0, new_line=.true.)
                write (6,'(1X,"#",1X,"Population on old reference det (averaged over report loop):",f10.2)') qs%D0_population
                write (6,'(1X,"#",1X,"Population on new reference det:",27X,i8)') max_pop
                write (6,'(1X,"#",1X,"E0 = <D0|H|D0> = ",f20.12)') reference%H00
                write (6,'(1X,"#",1X,"Care should be taken with accumulating statistics before this point.")')
                write (6,'(1X,"#",1X,62("-"))')
            end if
        end if

    end subroutine select_ref_det

    subroutine find_single_double_prob(sys, occ_list, psingle, pdouble)

        ! Calculate the probabilities of selecting a single or double excitation
        ! from a given determinant.  We assume all possible excitations (i.e.
        ! those with Hamiltonian matrix elements which are not zero by symmetry)
        ! are equally likely, so this amounts to finding the number of possible
        ! (symmetry-allowed) single and double excitations.
        !
        ! In:
        !    sys: system being studied.
        !    occ_list: integer list of occupied spin-orbitals in a determinant, D.
        ! Out:
        !    psingle: probability of attempting to spawn on a determinant
        !             connected to D by a single excitation.
        !    pdouble: probability of attempting to spawn on a determinant
        !             connected to D by a double excitation.

        use system
        use point_group_symmetry, only: cross_product_pg_basis, cross_product_pg_sym, nbasis_sym_spin

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: occ_list(sys%nel)
        real(p), intent(out) :: psingle, pdouble

        integer :: i, j, virt_syms(2, sys%sym0_tot:sys%sym_max_tot), nsingles, ndoubles, isyma, isymb, ims1, ims2

        select case(sys%system)
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
            do i = 1, sys%nel
                ! Convert -1->1 and 1->2 for spin index in arrays.
                ims1 = (sys%basis%basis_fns(occ_list(i))%ms+3)/2
                associate(isym=>sys%basis%basis_fns(occ_list(i))%sym)
                    virt_syms(ims1,isym) = virt_syms(ims1,isym) - 1
                end associate
            end do

            ! Count number of possible single excitations from the supplied
            ! determinant.
            ! Symmetry and spin must be conserved.
            nsingles = 0
            do i = 1, sys%nel
                ! Convert -1->1 and 1->2 for spin index in arrays.
                ims1 = (sys%basis%basis_fns(occ_list(i))%ms+3)/2
                ! Can't excite into already occupied orbitals.
                nsingles = nsingles + virt_syms(ims1,sys%basis%basis_fns(occ_list(i))%sym)
            end do

            ! Count number of possible double excitations from the supplied
            ! determinant.
            ndoubles = 0
            do i = 1, sys%nel
                ! Convert -1->1 and 1->2 for spin index in arrays.
                ims1 = (sys%basis%basis_fns(occ_list(i))%ms+3)/2
                do j = i+1, sys%nel
                    ! Convert -1->1 and 1->2 for spin index in arrays.
                    ims2 = (sys%basis%basis_fns(occ_list(j))%ms+3)/2
                    do isyma = sys%sym0, sys%sym_max
                        ! Symmetry of the final orbital is determined (for Abelian
                        ! symmetries) from the symmetry of the first three.
                        isymb = cross_product_pg_sym(isyma, cross_product_pg_basis(occ_list(i),occ_list(j), sys%basis%basis_fns))
                        if (isyma == isymb) then
                            if (ims1 == ims2) then
                                ! Cannot excit_t 2 electrons into the same spin-orbital.
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

    subroutine cumulative_population(pops, nactive, D0_proc, D0_pos, real_factor, cumulative_pops, tot_pop)

        ! Calculate the cumulative population, i.e. the number of psips/excips
        ! residing on a determinant/an excitor and all determinants/excitors which
        ! occur before it in the determinant/excitor list.

        ! This is primarily so in CCMC we can select clusters of excitors with each
        ! excip being equally likely to be part of a cluster.  (If we just select
        ! each occupied excitor with equal probability, then we get wildy
        ! fluctuating selection probabilities and hence large population blooms.)
        ! As 'excips' on the reference cannot be part of a cluster, then the
        ! population on the reference is treated as 0 if required.

        ! In:
        !    pops: list of populations on each determinant/excitor.  Must have
        !       minimum length of nactive.
        !    nactive: number of occupied determinants/excitors (ie pops(:,1:nactive)
        !       contains the population(s) on each currently "active"
        !       determinant/excitor.
        !    D0_proc: processor on which the reference resides.
        !    D0_pos: position in the pops list of the reference.  Only relevant if
        !       1<=D0_pos<=nactive and the processor holds the reference.
        !    real_factor: the encoding factor by which the stored populations are multiplied
        !       to enable non-integer populations.
        ! Out:
        !    cumulative_pops: running total of excitor population, i.e.
        !        cumulative_pops(i) = sum(abs(pops(1:i))), excluding the
        !        population on the reference if appropriate.
        !    tot_pop: total population (possibly excluding the population on the
        !       reference).

        ! NOTE: currently only the populations in the first psip/excip space are
        ! considered.  This should be changed if we do multiple simulations at
        ! once/Hellmann-Feynman sampling/etc.

        use parallel, only: iproc

        integer(int_p), intent(in) :: pops(:,:), real_factor
        integer, intent(in) :: nactive, D0_proc, D0_pos
        integer(int_p), intent(out) :: cumulative_pops(:), tot_pop

        integer :: i

! [review] - AJWT: This type of operation is almost certainly faster with a shift (or add and shift to get nint)
        cumulative_pops(1) = nint(real(abs(pops(1,1)),p)/real_factor)
        if (D0_proc == iproc) then
            ! Let's be a bit faster: unroll loops and skip over the reference
            ! between the loops.
            do i = 2, d0_pos-1
                cumulative_pops(i) = cumulative_pops(i-1) + &
                                        nint(real(abs(pops(1,i)),p)/real_factor)
            end do
            ! Set cumulative on the reference to be the running total merely so we
            ! can continue accessing the running total from the i-1 element in the
            ! loop over excitors in slots above the reference.
            if (d0_pos == 1) cumulative_pops(d0_pos) = 0
            if (d0_pos > 1) cumulative_pops(d0_pos) = cumulative_pops(d0_pos-1)
            do i = d0_pos+1, nactive
                cumulative_pops(i) = cumulative_pops(i-1) + &
                                        nint(real(abs(pops(1,i)),p)/real_factor)
            end do
        else
            ! V simple on other processors: no reference to get in the way!
            do i = 2, nactive
                cumulative_pops(i) = cumulative_pops(i-1) + &
                                        nint(real(abs(pops(1,i)),p)/real_factor)
            end do
        end if
        if (nactive > 0) then
            tot_pop = cumulative_pops(nactive)
        else
            tot_pop = 0
        end if

    end subroutine cumulative_population

    function decide_nattempts(rng, population) result(nattempts)

        ! Decide how many spawning attempts should be made from a determinant
        ! with the input population. If this population is not an integer, it
        ! must be stochastically rounded up or down in an unbiased manner.

        ! In:
        !    rng: random number generator.
        !    int_population: population of determinant, in its shifted integer
        !    form.

        use const, only: depsilon
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open

        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(in) :: population
        real(dp) :: r, pextra
        integer :: nattempts

        nattempts = abs(int(population))
        pextra = abs(population) - nattempts
        ! If there is no probability of generating an extra attempt, then
        ! don't bother using an extra random number.
        if (abs(pextra) > depsilon) then
            if (pextra > get_rand_close_open(rng)) nattempts = nattempts + 1
        end if

    end function decide_nattempts

    subroutine load_balancing_report(nparticles, nstates_active, spawn_mpi_time, determ_mpi_time)

        ! In:
        !    nparticles: number of particles in each space, on this process only.
        !    nstates_active: number of occupied states, on this process only.
        !    spawn_mpi_time: MPI timings for the spawned list, on this process
        !        only.
        ! In (optional):
        !    determ_mpi_time: MPI timings for semi-stochastic communications,
        !        on this process only.

        ! Print out a load-balancing report when run in parallel showing how
        ! determinants and walkers/particles are distributed over the processors.

#ifdef PARALLEL
        use calc, only: use_mpi_barriers
        use parallel
        use utils, only: int_fmt

        real(p), intent(in) :: nparticles(:)
        integer, intent(in) :: nstates_active
        type(parallel_timing_t), intent(in) :: spawn_mpi_time
        type(parallel_timing_t), optional, intent(in) :: determ_mpi_time

        real(dp) :: load_data(size(nparticles), nprocs)
        integer :: load_data_int(nprocs)
        integer :: i, ierr
        real(p) :: barrier_this_proc
        real(p) :: spawn_comms(nprocs), determ_comms(nprocs), barrier_time(nprocs)
        character(4) :: lfmt

        if (nprocs > 1) then
            if (parent) then
                write (6,'(1X,a14,/,1X,14("^"),/)') 'Load balancing'
                write (6,'(1X,a77,/)') "The final distribution of walkers and determinants across the processors was:"
            endif
            call mpi_gather(nparticles, size(nparticles), mpi_preal, load_data, size(nparticles), &
                            mpi_preal, 0, MPI_COMM_WORLD, ierr)
            if (parent) then
                do i = 1, size(nparticles)
                    if (size(nparticles) > 1) write (6,'(1X,a,'//int_fmt(i,1)//')') 'Particle type:', i
                    write (6,'(1X,"Min # of particles on a processor:",6X,es12.6)') minval(load_data(i,:))
                    write (6,'(1X,"Max # of particles on a processor:",6X,es12.6)') maxval(load_data(i,:))
                    write (6,'(1X,"Mean # of particles on a processor:",5X,es12.6,/)') real(sum(load_data(i,:)), p)/nprocs
                end do
            end if
            call mpi_gather(nstates_active, 1, mpi_integer, load_data_int, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr)
            call mpi_gather(spawn_mpi_time%comm_time, 1, mpi_preal, spawn_comms, 1, mpi_preal, 0, MPI_COMM_WORLD, ierr)

            if (present(determ_mpi_time)) call mpi_gather(determ_mpi_time%comm_time, 1, mpi_preal, determ_comms, 1, &
                                                           mpi_preal, 0, MPI_COMM_WORLD, ierr)

            if (use_mpi_barriers) then
                if (present(determ_mpi_time)) then
                    barrier_this_proc = spawn_mpi_time%barrier_time + determ_mpi_time%barrier_time
                else
                    barrier_this_proc = spawn_mpi_time%barrier_time
                end if
                call mpi_gather(barrier_this_proc, 1, mpi_preal, barrier_time, 1, mpi_preal, 0, MPI_COMM_WORLD, ierr)
            end if

            if (parent) then
                lfmt = int_fmt(maxval(load_data_int),0)
                write (6,'(1X,"Min # of determinants on a processor:",3X,'//lfmt//')') minval(load_data_int)
                write (6,'(1X,"Max # of determinants on a processor:",3X,'//lfmt//')') maxval(load_data_int)
                write (6,'(1X,"Mean # of determinants on a processor:",2X,es12.6)') real(sum(load_data_int), p)/nprocs
                write (6,'()')
                if (use_mpi_barriers) then
                    write (6,'(1X,"Min time taken by MPI barrier calls:",5X,f8.2,"s")') minval(barrier_time)
                    write (6,'(1X,"Max time taken by MPI barrier calls:",5X,f8.2,"s")') maxval(barrier_time)
                    write (6,'(1X,"Mean time taken by MPI barrier calls:",4X,f8.2,"s")') real(sum(barrier_time), p)/nprocs
                    write (6,'()')
                end if
                write (6,'(1X,"Min time taken by walker communication:",5X,f8.2,"s")') minval(spawn_comms)
                write (6,'(1X,"Max time taken by walker communication:",5X,f8.2,"s")') maxval(spawn_comms)
                write (6,'(1X,"Mean time taken by walker communication:",4X,f8.2,"s")') real(sum(spawn_comms), p)/nprocs
                write (6,'()')
                if (present(determ_mpi_time)) then
                    write (6,'(1X,"Min time taken by semi-stochastic communication:",5X,f8.2,"s")') minval(determ_comms)
                    write (6,'(1X,"Max time taken by semi-stochastic communication:",5X,f8.2,"s")') maxval(determ_comms)
                    write (6,'(1X,"Mean time taken by semi-stochastic communication:",4X,f8.2,"s")') &
                        real(sum(determ_comms), p)/nprocs
                    write (6,'()')
                end if
            end if
        end if
#else
        use parallel, only: parallel_timing_t

        type(parallel_timing_t), intent(in) :: spawn_mpi_time
        type(parallel_timing_t), optional, intent(in) :: determ_mpi_time
#endif

    end subroutine load_balancing_report

    subroutine redistribute_particles(states, real_factor, pops, nstates, nparticles, spawn, &
                                      nload_slots)

        ! [todo] JSS: - update comments to be more general than just for CCMC.

        ! Due to the cooperative spawning (ie from multiple excitors at once) in
        ! CCMC, we need to give each excitor the chance to be on the same
        ! processor with all combinations of excitors, unlike in FCIQMC where
        ! the spawning events are independent.  We satisfy this by periodically
        ! moving an excitor to a different processor (MPI rank).

        ! WARNING: if the number of processors is large or the system small,
        ! this introduces a bias as load balancing prevents all possible
        ! clusters from being on the same processor at the same time.

        ! In:
        !    states: list of occupied excitors on the current processor.
        !    real_factor: The factor by which populations are multiplied to
        !        enable non-integer populations.
        !    nload_slots: number of load balancing slots (per processor).
        ! In/Out:
        !    nparticles: number of excips on the current processor.
        !    pops: Population on occupied excitors.  On output the
        !        populations of excitors which are sent to other processors are
        !        set to zero.
        !    nstates: number of occupied excitors on the current processor.
        !    spawn: spawn_t object.  On output particles which need to be sent
        !        to another processor have been added to the correct position in
        !        the spawned store.

        use const, only: i0, dp
        use spawn_data, only: spawn_t
        use spawning, only: assign_particle_processor, add_spawned_particles
        use parallel, only: iproc, nprocs

        integer(i0), intent(in) :: states(:,:)
        integer(int_p), intent(in) :: real_factor
        integer(int_p), intent(inout) :: pops(:,:)
        integer, intent(inout) :: nstates
        real(p), intent(inout) :: nparticles(:)
        type(spawn_t), intent(inout) :: spawn
        integer, intent(in) :: nload_slots

        real(p) :: nsent(size(nparticles))

        integer :: iexcitor, pproc, string_len, slot

        nsent = 0.0_dp
        string_len = size(states, dim=1)

        !$omp parallel do default(none) &
        !$omp shared(nstates, %states, %popss, spawn, iproc, nprocs, string_len) &
        !$omp private(pproc, slot) reduction(+:nsent)
        do iexcitor = 1, nstates
            !  - set hash_shift and move_freq
            call assign_particle_processor(states(:,iexcitor), string_len, spawn%hash_seed, &
                                           spawn%hash_shift, spawn%move_freq, nprocs, pproc, slot, nload_slots)
            if (pproc /= iproc) then
                ! Need to move.
                ! Add to spawned array so it will be sent to the correct
                ! processor during annihilation.
                ! NOTE: for initiator calculations we need to keep this
                ! population no matter what.  This relies upon the
                ! (undocumented) 'feature' that a flag of 0 indicates the parent
                ! was an initiator...
                call add_spawned_particles(states(:,iexcitor), pops(:,iexcitor), pproc, spawn)
                ! Update population on the sending processor.
                nsent = nsent + abs(real(pops(:,iexcitor),p))
                ! Zero population here.  Will be pruned on this determinant
                ! automatically during annihilation (which will also update nstates).
                pops(:,iexcitor) = 0_int_p
            end if
        end do
        !$omp end parallel do

        ! Remove encoding factor to obtain the true populations.
        nsent = nsent/real_factor

        nparticles = nparticles - nsent

    end subroutine redistribute_particles

    subroutine redistribute_semi_stoch_t(sys, reference, annihilation_flags, psip_list, spawn, determ, nload_slots)

        ! Recreate the semi_stoch_t object (if a non-empty space is in use).
        ! This requires sending deterministic states to their new processes
        ! and recreating the related objects, such as the deterministic
        ! Hamiltonian.

        ! In:
        !    sys: system being studied.
        !    annihilation_flags: calculation specific annihilation flags.
        !    psip_list: list of particles and their locations.
        !    spawn: spawn_t object, required for determining the new processes
        !        labels for deterministic states.
        !    nload_slots: number of load balancing slots (per processor).
        ! In/Out:
        !    determ: The deterministic space being used, as required for
        !        semi-stochastic calculations.

        use semi_stoch, only: dealloc_semi_stoch_t, init_semi_stoch_t

        use spawn_data, only: spawn_t
        use system, only: sys_t
        use qmc_data, only: reference_t, semi_stoch_t, empty_determ_space, reuse_determ_space, particle_t
        use qmc_data, only: annihilation_flags_t

        type(sys_t), intent(in) :: sys
        type(annihilation_flags_t), intent(in) :: annihilation_flags
        type(reference_t), intent(in) :: reference
        type(particle_t), intent(inout) :: psip_list
        type(spawn_t), intent(in) :: spawn
        type(semi_stoch_t), intent(inout) :: determ
        integer, intent(in) :: nload_slots

        logical :: sep_annihil_copy

        if (determ%space_type /= empty_determ_space) then
            ! Copy this logical to re-pass into the initialisation routine.
            sep_annihil_copy = determ%separate_annihilation

            ! Deallocate the semi_stoch_t instance, except for the list of all
            ! deterministic states, which we want to reuse.
            call dealloc_semi_stoch_t(determ, keep_dets=.true.)
            ! Recreate the semi_stoch_t instance, by reusing the deterministic
            ! space already generated, but with states on their new processes.
            call init_semi_stoch_t(determ, sys, psip_list, reference, annihilation_flags, spawn, reuse_determ_space, &
                                    0, sep_annihil_copy, .false., .true., nload_slots)
        end if

    end subroutine redistribute_semi_stoch_t

! --- Output routines ---

    subroutine initial_fciqmc_status(sys, qmc_in, qs, reference, psip_list, nb_comm, rep_comm, spawn_elsewhere)

        ! Calculate the projected energy based upon the initial walker
        ! distribution (either via a restart or as set during initialisation)
        ! and print out.

        ! In:
        !    sys: system being studied.
        !    qmc_in: input options relating to QMC methods.
        !    reference: current reference determinant.
        !    psip_list: set of particles representing the current representation
        !       of the section of the wavefunction stored  on the processor.
        ! In (optional):
        !    nb_comm: true if using non-blocking communications.
        !    spawn_elsewhere: number of particles spawned from the current
        !       processor to other processors.  Relevant only when restarting
        !       non-blocking calculations.
        ! Out (Optional):
        !    rep_comm: nb_rep_t object containg report loop information.

        use calc, only: nb_rep_t
        use determinants, only: det_info_t, alloc_det_info_t, dealloc_det_info_t, decode_det
        use energy_evaluation, only: local_energy_estimators, update_energy_estimators_send
        use excitations, only: excit_t, get_excitation
        use parallel
        use proc_pointers, only: update_proj_energy_ptr
        use qmc_data, only: qmc_in_t, reference_t, particle_t, qmc_state_t
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        type(qmc_state_t), intent(inout) :: qs
        type(reference_t), intent(in) :: reference
        type(particle_t), intent(in), target :: psip_list
        logical, optional, intent(in) :: nb_comm
        type(nb_rep_t), optional, intent(inout) :: rep_comm
        integer, optional, intent(in) :: spawn_elsewhere

        integer :: idet
        real(p) :: ntot_particles(psip_list%nspaces)
        real(p) :: real_population(psip_list%nspaces)
        type(det_info_t) :: cdet
        real(p) :: hmatel
        type(excit_t) :: D0_excit
        logical :: nb_comm_local
#ifdef PARALLEL
        integer :: ierr
        real(p) :: proj_energy_sum, D0_population_sum
#endif

        ! Calculate the projected energy based upon the initial walker
        ! distribution.  proj_energy and D0_population are both accumulated in
        ! update_proj_energy.
        qs%proj_energy = 0.0_p
        qs%D0_population = 0.0_p
        call alloc_det_info_t(sys, cdet)
        do idet = 1, psip_list%nstates
            cdet%f = psip_list%states(:,idet)
            call decode_det(sys%basis, cdet%f, cdet%occ_list)
            cdet%data => psip_list%dat(:,idet)
            real_population = real(psip_list%pops(:,idet),p)/real_factor
            ! WARNING!  We assume only the bit string, occ list and data field
            ! are required to update the projected estimator.
            D0_excit = get_excitation(sys%nel, sys%basis, cdet%f, reference%f0)
            call update_proj_energy_ptr(sys, reference%f0, cdet, real_population(1), &
                                        qs%D0_population, qs%proj_energy, D0_excit, hmatel)
        end do
        call dealloc_det_info_t(cdet)

#ifdef PARALLEL
        if (present(rep_comm)) then
            ! The output in non-blocking comms is delayed one report loop, so initialise
            ! the send here.
            ! For simplicity, hook into the normal estimator communications, which normalises
            ! by the number of MC cycles in a report loop (hence need to rescale to fake it).
            qs%D0_population = qs%D0_population*qmc_in%ncycles
            qs%proj_energy = qs%proj_energy*qmc_in%ncycles
            call local_energy_estimators(psip_list, qs, rep_comm%rep_info, spawn_elsewhere=spawn_elsewhere)
            call update_energy_estimators_send(rep_comm)
        else
            call mpi_allreduce(qs%proj_energy, proj_energy_sum, psip_list%nspaces, mpi_preal, &
                               MPI_SUM, MPI_COMM_WORLD, ierr)
            call mpi_allreduce(psip_list%nparticles, ntot_particles, psip_list%nspaces, MPI_PREAL, &
                               MPI_SUM, MPI_COMM_WORLD, ierr)
            call mpi_allreduce(qs%D0_population, D0_population_sum, 1, mpi_preal, MPI_SUM, MPI_COMM_WORLD, ierr)
            qs%proj_energy = proj_energy_sum
            qs%D0_population = D0_population_sum
            ! TODO: HFS, DMQMC quantities
        end if
#else
        ntot_particles = psip_list%nparticles
#endif

        ! Using non blocking communications?
        nb_comm_local = .false.
        if (present(nb_comm)) nb_comm_local = nb_comm

        if (.not. nb_comm_local .and. parent) then
            ! See also the format used in write_fciqmc_report if this is changed.
            ! We prepend a # to make it easy to skip this point when do data
            ! analysis.
            call write_fciqmc_report(qmc_in, qs, 0, ntot_particles, 0.0, .true., .false.)
        end if

    end subroutine initial_fciqmc_status

! --- QMC loop and cycle initialisation routines ---

    subroutine init_report_loop(qs, bloom_stats)

        ! Initialise a report loop (basically zero quantities accumulated over
        ! a report loop).

        use bloom_handler, only: bloom_stats_t, bloom_stats_init_report_loop
        use qmc_data, only: qmc_state_t

        type(qmc_state_t), intent(inout) :: qs
        type(bloom_stats_t), intent(inout) :: bloom_stats

        call bloom_stats_init_report_loop(bloom_stats)

        qs%proj_energy = 0.0_p
        rspawn = 0.0_p
        qs%D0_population = 0.0_p

        ! DMQMC-specific...
        if (allocated(excit_dist)) excit_dist = 0.0_p
        if (allocated(trace)) trace = 0.0_p
        numerators = 0.0_p

    end subroutine init_report_loop

    subroutine init_mc_cycle(rng, sys, qmc_in, reference, load_bal_in, annihilation_flags, real_factor, &
                             psip_list, spawn, nattempts, ndeath, min_attempts, doing_lb, nb_comm, determ)

        ! Initialise a Monte Carlo cycle (basically zero/reset cycle-level
        ! quantities).

        ! In/Out:
        !    rng: random number generator.
        ! In:
        !    sys: system being studied
        !    qmc_in: input options relating to QMC methods.
        !    reference: current reference determinant.
        !    load_bal_in: input options for load balancing.
        !    annihilation_flags: calculation specific annihilation flags.
        !    real_factor: The factor by which populations are multiplied to
        !        enable non-integer populations.
        ! In/Out:
        !    psip_list: total population (on this proccesor) is used to set
        !       nattempts and population is redistributed if requested by the
        !       load balancing approach.
        !    spawn: spawn_t object for holding spawned particles.  Reset on exit.
        ! Out:
        !    nattempts: number of spawning attempts to be made (on the current
        !        processor) this cycle.
        !    ndeath: number of particle deaths that occur in a Monte Carlo
        !        cycle.  Reset to 0 on output.
        ! In (optional):
        !    min_attempts: if present, set nattempts to be at least this value.
        !    doing_lb: true if doing load balancing.
        !    nb_comm: true if using non-blocking communications.
        !    nload_slots: number of load balancing slots proc_map is divided
        !        into.
        ! In/Out (optional):
        !    determ: The deterministic space being used, as required for
        !        semi-stochastic calculations.

        use calc, only: doing_calc, ct_fciqmc_calc, ccmc_calc, dmqmc_calc
        use dSFMT_interface, only: dSFMT_t
        use load_balancing, only: do_load_balancing
        use qmc_data, only: qmc_in_t, reference_t, load_bal_in_t, semi_stoch_t, particle_t, annihilation_flags_t
        use system, only: sys_t
        use spawn_data, only: spawn_t

        type(dSFMT_t), intent(inout) :: rng
        type(sys_t), intent(in) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        type(reference_t), intent(in) :: reference
        type(annihilation_flags_t), intent(in) :: annihilation_flags
        type(particle_t), intent(inout) :: psip_list
        type(spawn_t), intent(inout) :: spawn
        type(load_bal_in_t), intent(in) :: load_bal_in
        integer(int_p), intent(in) :: real_factor
        integer(int_64), intent(in), optional :: min_attempts
        integer(int_64), intent(out) :: nattempts
        integer(int_p), intent(out) :: ndeath
        logical, optional, intent(in) :: doing_lb, nb_comm
        type(semi_stoch_t), optional, intent(inout) :: determ

        logical :: nb_comm_local

        ! Using non-blocking communications?
        nb_comm_local = .false.
        if (present(nb_comm)) nb_comm_local = nb_comm

        ! Reset the current position in the spawning array to be the
        ! slot preceding the first slot.
        spawn%head = spawn%head_start

        ! Reset death counter
        ndeath = 0_int_p

        ! Number of spawning attempts that will be made.
        ! For FCIQMC, this is used for accounting later, not for controlling the
        ! spawning.
        if (doing_calc(ct_fciqmc_calc) .or. doing_calc(ccmc_calc)) then
            ! ct algorithm: kinda poorly defined.
            ! ccmc: number of excitor clusters we'll randomly generate and
            ! attempt to spawn from.
            nattempts = psip_list%nparticles(1)
        else if (doing_calc(dmqmc_calc)) then
            ! Each particle and each end gets to attempt to spawn onto a
            ! connected determinant and a chance to die/clone.
            nattempts = nint(4*psip_list%nparticles(1)*psip_list%nspaces)
        else
            ! Each particle gets to attempt to spawn onto a connected
            ! determinant and a chance to die/clone.
            nattempts = nint(2*psip_list%nparticles(1))
        end if

        if (present(min_attempts)) nattempts = max(nattempts, min_attempts)

        if (present(doing_lb)) then
            if (doing_lb .and. par_info%load%needed) then
                call do_load_balancing(psip_list, spawn, real_factor, par_info, load_bal_in)
                call redistribute_load_balancing_dets(rng, sys, qmc_in, reference, psip_list%states, real_factor, determ, &
                                                      psip_list, spawn, load_bal_in%nslots, annihilation_flags)
                ! If using non-blocking communications we still need this flag to
                ! be set.
                if (.not. nb_comm_local) par_info%load%needed = .false.
            end if
        end if

    end subroutine init_mc_cycle

! --- QMC loop and cycle termination routines ---

    subroutine end_report_loop(sys, qmc_in, reference, ireport, iteration, update_tau, psip_list, qs,  &
                                ntot_particles, nspawn_events, report_time, semi_stoch_shift_it,       &
                                semi_stoch_start_it, soft_exit, dump_restart_file_shift, load_bal_in,  &
                                update_estimators, bloom_stats, doing_lb, nb_comm, rep_comm, dmqmc_in, &
                                spawn_recv)

        ! In:
        !    sys: system being studied.
        !    reference: current reference determinant.
        !    ireport: index of current report loop.
        !    iteration: The current iteration of the simulation.
        !    update_tau: true if the processor thinks the timestep should be rescaled.
        !             Only used if not in variable shift mode and if tau_search is being
        !             used.
        !    nspawn_events: The total number of spawning events to this process.
        !    semi_stoch_shift_it: How many iterations after the shift starts
        !        to vary to begin using semi-stochastic.
        !    load_bal_in: input options for load balancing.
        !    
        ! In/Out:
        !    qmc_in: input optons relating to QMC methods.
        !    psip_list: particle_t object with current psip distribution on the
        !        current proccesor.  Total (global) quantities are updated on
        !        exit.
        !    ntot_particles: total number (across all processors) of
        !        particles in the simulation at end of the previous report loop.
        !        Returns the current total number of particles for use in the
        !        next report loop if update_estimators is true.
        !    semi_stoch_start_it: The iteration on which to start performing
        !        semi-stochastic.
        !    report_time: time at the start of the current report loop.  Returns
        !        the current time (ie the time for the start of the next report
        !        loop.
        !    dump_restart_file_shift: should we dump a restart file just before
        !        the shift turns on?
        ! Out:
        !    soft_exit: true if the user has requested an immediate exit of the
        !        QMC algorithm via the interactive functionality.
        ! In (optional):
        !    update_estimators: update the (FCIQMC/CCMC) energy estimators.  Default: true.
        !    doing_lb: true if doing load balancing.
        !    nb_comm: true if using non-blocking communications.
        !    load_bal_in: input options for load balancing.
        !    dmqmc_in: input options relating to DMQMC.
        !    spawn_recv: spawn_t object containing particles spawned in the previous iteration and received using non-blocking comms
        !                algorithm.
        ! In/Out (optional):
        !    bloom_stats: particle blooming statistics to accumulate.
        !    rep_comm: nb_rep_t object containing report loop info. Used for
        !        non-blocking communications where we receive report information
        !        from previous iteration and communicate the current iterations
        !        estimators.

        use dmqmc_data, only: dmqmc_in_t
        use energy_evaluation, only: update_energy_estimators, local_energy_estimators,         &
                                     update_energy_estimators_recv, update_energy_estimators_send, &
                                     nparticles_start_ind
        use interact, only: calc_interact, check_interact, check_comms_file
        use parallel, only: parent, nprocs
        use restart_hdf5, only: dump_restart_hdf5, restart_info_global, restart_info_global_shift
        use system, only: sys_t
        use calc, only: nb_rep_t
        use bloom_handler, only: bloom_stats_t, bloom_stats_warning
        use qmc_data, only: qmc_in_t, reference_t, load_bal_in_t, particle_t, qmc_state_t

        type(sys_t), intent(in) :: sys
        type(reference_t), intent(in) :: reference
        type(qmc_in_t), intent(inout) :: qmc_in
        integer, intent(in) :: ireport, iteration
        logical, intent(in) :: update_tau
        type(particle_t), intent(inout) :: psip_list
        type(qmc_state_t), intent(inout) :: qs
        integer, intent(in) :: nspawn_events
        logical, optional, intent(in) :: update_estimators
        type(bloom_stats_t), optional, intent(inout) :: bloom_stats
        real(p), intent(inout) :: ntot_particles(psip_list%nspaces)
        real, intent(inout) :: report_time
        integer, intent(in) :: semi_stoch_shift_it
        integer, intent(inout) :: semi_stoch_start_it
        logical, intent(out) :: soft_exit
        logical, intent(inout) :: dump_restart_file_shift
        type(load_bal_in_t), intent(in) :: load_bal_in
        logical, optional, intent(in) :: doing_lb, nb_comm
        type(nb_rep_t), optional, intent(inout) :: rep_comm
        type(dmqmc_in_t), optional, intent(in) :: dmqmc_in
        type(spawn_t), optional, intent(in) :: spawn_recv

        real :: curr_time
        logical :: update, update_tau_now, vary_shift_before, nb_comm_local, comms_found
        real(dp) :: rep_info_copy(nprocs*psip_list%nspaces+nparticles_start_ind-1)

        ! Only update the timestep if not in vary shift mode.
        update_tau_now = update_tau .and. .not. qs%vary_shift(1) .and. qmc_in%tau_search

        ! Using non-blocking communications?
        nb_comm_local = .false.
        if (present(nb_comm)) nb_comm_local = nb_comm

        ! Are all the shifts currently varying?
        vary_shift_before = all(qs%vary_shift)

        ! Test for a comms file so MPI communication can be combined with
        ! energy_estimators communication
        comms_found = check_comms_file()

        ! Update the energy estimators (shift & projected energy).
        update = .true.
        if (present(update_estimators)) update = update_estimators
        if (update .and. .not. nb_comm_local) then
            call update_energy_estimators(qmc_in, qs, nspawn_events, psip_list, ntot_particles, load_bal_in, doing_lb, comms_found, &
                                          update_tau_now, bloom_stats)
        else if (update) then
            ! Save current report loop quantitites.
            ! Can't overwrite the send buffer before message completion
            ! so copy information somewhere else.
            call local_energy_estimators(psip_list, qs, rep_info_copy, nspawn_events, comms_found, update_tau_now, bloom_stats, &
                                          rep_comm%nb_spawn(2))
            ! Receive previous iterations report loop quantities.
            call update_energy_estimators_recv(qmc_in, qs, psip_list%nspaces, rep_comm%request, ntot_particles,             &
                                               psip_list%nparticles_proc, load_bal_in, doing_lb, comms_found, update_tau_now, &
                                               bloom_stats)
            ! Send current report loop quantities.
            rep_comm%rep_info = rep_info_copy
            call update_energy_estimators_send(rep_comm)
        else
            update_tau_now = .false.
            call check_interact(comms_found)
        end if

        ! If we have just started varying the shift, then we can calculate the
        ! iteration at which to start semi-stochastic, if it is being turned on
        ! relative to the shift start.
        if ((.not. vary_shift_before) .and. all(qs%vary_shift) .and. (semi_stoch_shift_it /= -1)) &
            semi_stoch_start_it = semi_stoch_shift_it + iteration + 1

        call cpu_time(curr_time)

        ! report_time was the time at the previous iteration.
        ! curr_time - report_time is thus the time taken by this report loop.
        if (parent) then
            if (bloom_stats%nblooms_curr > 0) call bloom_stats_warning(bloom_stats)
            call write_fciqmc_report(qmc_in, qs, ireport, ntot_particles, curr_time-report_time, .false., nb_comm_local, dmqmc_in)
        end if

        ! Write restart file if required.
        if (dump_restart_file_shift .and. any(qs%vary_shift)) then
            dump_restart_file_shift = .false.
            call dump_restart_hdf5(restart_info_global_shift, qs, psip_list, reference, mc_cycles_done+qmc_in%ncycles*ireport, &
                                   ntot_particles, nb_comm_local, spawn_recv)
        else if (mod(ireport,restart_info_global%write_restart_freq) == 0) then
            call dump_restart_hdf5(restart_info_global, qs, psip_list, reference, mc_cycles_done+qmc_in%ncycles*ireport, &
                                   ntot_particles, nb_comm_local, spawn_recv)
        end if
        ! cpu_time outputs an elapsed time, so update the reference timer.
        report_time = curr_time

        call calc_interact(comms_found, soft_exit, qmc_in)

        if (update_tau_now) call rescale_tau(qs%tau)

    end subroutine end_report_loop

    subroutine end_mc_cycle(nspawn_events, ndeath, nattempts)

        ! Execute common code at the end of a Monte Carlo cycle.

        ! In:
        !    nspawn_events: number of successful spawning events in the cycle.
        !    ndeath: (unscaled) number of particle deaths in the cycle.
        !    nattempts: number of attempted spawning events in the cycle.

        integer, intent(in) :: nspawn_events
        integer(int_p), intent(in) :: ndeath
        integer(int_64), intent(in) :: nattempts

        ! Add the spawning rate (for the processor) to the running
        ! total.
        rspawn = rspawn + spawning_rate(nspawn_events, ndeath, nattempts)

    end subroutine end_mc_cycle

    subroutine rescale_tau(tau, factor)

        ! Scale the timestep by across all the processors.  This is performed if at least
        ! one processor thinks it should be.

        ! In:
        !    tau: timestep to be updated.
        !    factor: factor to scale the timestep by (default: 0.95).

        use parallel

        real(p), intent(inout) :: tau
        real(p), intent(in), optional :: factor

        integer :: ierr
        logical :: update_tau_global

        if (present(factor)) then
            tau = factor*tau
        else
            tau = 0.950_p*tau
        end if
        if (parent) write(6, '(1X, "# Warning timestep changed to: ",f7.5)') tau

    end subroutine rescale_tau

    subroutine redistribute_load_balancing_dets(rng, sys, qmc_in, reference, states, real_factor, &
                                                determ, psip_list, spawn, nload_slots, annihilation_flags)

        ! When doing load balancing we need to redistribute chosen sections of
        ! main list to be sent to their new processors. This is a wrapper which
        ! takes care of this and resets the load balancing tag so that no
        ! further load balancing is attempted this report loop. Currently don't
        ! check if anything will actually be sent.

        ! Also if a non-empty semi-stochastic deterministic space is being used
        ! then this object needs to be redistributed, too.

        ! In:
        !    sys: system being studied.
        !    qmc_in: input options relating to QMC methods.
        !    reference: current reference determinant.
        !    states: list of occupied excitors on the current processor.
        !    real_factor: The factor by which populations are multiplied to
        !        enable non-integer populations.
        !    nload_slots: number of load balancing slots (per processor).
        !    annihilation_flags: calculation specific annihilation flags.
        ! In/Out:
        !    rng: random number generator.
        !    determ (optional): The deterministic space being used, as required for
        !        semi-stochastic calculations.
        !    psip_list: particle_t object containing current distribution of
        !       psips.  On exit the list contains a set which (in principle)
        !       improves load balancing by sending/receiving particles from
        !       other processors.  All components are updated as required for
        !       consistency.
        !    spawn: spawn_t object.  Used to send particles to their new
        !       processor.

        use annihilation, only: direct_annihilation
        use dSFMT_interface, only: dSFMT_t
        use qmc_data, only: qmc_in_t, reference_t, semi_stoch_t, particle_t, annihilation_flags_t
        use spawn_data, only: spawn_t
        use system, only: sys_t

        type(dSFMT_t), intent(inout) :: rng
        type(sys_t), intent(in) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        type(reference_t), intent(in) :: reference
        integer(i0), intent(in) :: states(:,:)
        integer(int_p), intent(in) :: real_factor
        type(semi_stoch_t), optional, intent(inout) :: determ
        type(particle_t), intent(inout) :: psip_list
        type(spawn_t), intent(inout) :: spawn
        integer, intent(in) :: nload_slots
        type(annihilation_flags_t), intent(in) :: annihilation_flags

        call redistribute_particles(psip_list%states, real_factor, psip_list%pops, &
                                    psip_list%nstates, psip_list%nparticles, spawn, nload_slots)

        ! Merge determinants which have potentially moved processor back into
        ! the appropriate main list.
        call direct_annihilation(sys, rng, qmc_in, reference, annihilation_flags, psip_list, spawn)
        spawn%head = spawn%head_start

        if (present(determ)) call redistribute_semi_stoch_t(sys, reference, annihilation_flags, psip_list, &
                                                            spawn, determ, nload_slots)

    end subroutine redistribute_load_balancing_dets

end module qmc_common
