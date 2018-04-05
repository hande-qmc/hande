module qmc_common

! Module containing routines common to different qmc algorithms.

use const

implicit none

contains

! --- Utility routines ---

    subroutine select_ref_det(sys, ref_det_factor, qs)

        ! Change the reference determinant to be the determinant with the
        ! greatest population if it exceeds some threshold relative to the
        ! current reference determinant.

        ! Note this looks at the sum of populations over all spaces.  The
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
        !    qs: qmc_state_t object. On output the reference determinant may be
        !        updated, and the diagagonal matrix elements in psip_list also.

        use calc, only: doing_calc, hfs_fciqmc_calc
        use determinants, only: decode_det, write_det, sum_sp_eigenvalues_occ_list
        use system, only: sys_t

        use parallel
        use errors, only: stop_all
        use qmc_data, only: qmc_state_t

        type(sys_t), intent(in) :: sys
        real(p), intent(in) :: ref_det_factor
        type(qmc_state_t), intent(inout) :: qs

        integer :: i
        integer(i0), allocatable :: fmax(:)
        integer(int_p) :: max_pop
#ifdef PARALLEL
        real(dp) :: in_data(2), out_data(2)
        integer :: D0_proc, ierr
#endif
        real(p) :: H00_max, H00_old
        real(dp) :: real_pop
        logical :: updated
        integer :: iunit

        iunit = 6

        allocate(fmax(lbound(qs%psip_list%states, dim=1):ubound(qs%psip_list%states, dim=1)))

        H00_old = qs%ref%H00

        updated = .false.
        ! Find determinant with largest population.
        max_pop = 0_int_p
        do i = 1, qs%psip_list%nstates
            if (sum(abs(qs%psip_list%pops(:,i))) > abs(max_pop)) then
                max_pop = sum(abs(qs%psip_list%pops(:,i)))
                fmax = qs%psip_list%states(:,i)
                H00_max = qs%psip_list%dat(1, i)
            end if
        end do

        real_pop = real(max_pop,dp)/qs%psip_list%pop_real_factor

        ! Only change reference determinant if the population is larger than the
        ! reference determinant by a given factor to avoid switching
        ! continuously between degenerate determinants.

        ! Note we don't broadcast the population of the new reference det as
        ! that is reset at the start of the next report loop anyway (and this
        ! routine should only be called at the end of the report loop).

#ifdef PARALLEL

        if (all(fmax == qs%ref%f0)) then
            ! Max population on this processor is already the reference.  Don't change.
            in_data = (/ 0.0_dp, real(iproc,dp) /)
        else if (abs(real_pop) > ref_det_factor*sum(abs(qs%estimators%D0_population))) then
            in_data = (/ real_pop, real(iproc,dp) /)
        else
            ! No det with sufficient population to become reference det on this
            ! processor.
            in_data = (/ 0.0_dp, real(iproc, dp) /)
        end if

        call mpi_allreduce(in_data, out_data, 1, mpi_2double_precision, MPI_MAXLOC, MPI_COMM_WORLD, ierr)

        if (abs(out_data(1)) > depsilon) then
            real_pop = out_data(1)
            updated = .true.
            D0_proc = nint(out_data(2))
            qs%ref%f0 = fmax
            qs%ref%H00 = H00_max
            ! Broadcast updated data
            call mpi_bcast(qs%ref%f0, size(qs%ref%f0), mpi_det_integer, D0_proc, MPI_COMM_WORLD, ierr)
            call mpi_bcast(qs%ref%H00, 1, mpi_preal, D0_proc, MPI_COMM_WORLD, ierr)
        end if

#else

        if (abs(real_pop) > ref_det_factor*sum(abs(qs%estimators%D0_population)) .and. any(fmax /= qs%ref%f0)) then
            updated = .true.
            qs%ref%f0 = fmax
            qs%ref%H00 = H00_max
        end if

#endif

        if (updated) then
            ! Update occ_list.
            call decode_det(sys%basis, qs%ref%f0, qs%ref%occ_list0)
            ! psip_list%dat(1,i) holds <D_i|H|D_i> - H00_old.  Update.
            ! H00 is currently <D_0|H|D_0> - H00_old.
            ! Want psip_list%dat(1,i) to be <D_i|H|D_i> - <D_0|H|D_0>
            ! We'll fix H00 later and avoid an extra nstates*additions.
            do i = 1, qs%psip_list%nstates
                qs%psip_list%dat(1,i) = qs%psip_list%dat(1,i) - qs%ref%H00
            end do
            ! Now set H00 = <D_0|H|D_0> so that future references to it are
            ! correct.
            qs%ref%H00 = qs%ref%H00 + H00_old
            qs%ref%fock_sum = sum_sp_eigenvalues_occ_list(sys, qs%ref%occ_list0)
            if (doing_calc(hfs_fciqmc_calc)) call stop_all('select_ref_det', 'Not implemented for HFS.')
            if (parent) then
                write (iunit,'(1X,"#",1X,62("-"))')
                write (iunit,'(1X,"#",1X,"Changed reference det to:",1X)',advance='no')
                call write_det(sys%basis, sys%nel, qs%ref%f0, new_line=.true.)
                write (iunit,'(1X,"#",1X,"Population on old reference det (averaged over report loop):",f10.2)') &
                            sum(abs(qs%estimators%D0_population))
                write (iunit,'(1X,"#",1X,"Population on new reference det:",27X,f10.2)') real_pop
                write (iunit,'(1X,"#",1X,"E0 = <D0|H|D0> = ",f20.12)') qs%ref%H00
                write (iunit,'(1X,"#",1X,"Care should be taken with accumulating statistics before this point.")')
                write (iunit,'(1X,"#",1X,62("-"))')
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
        use symmetry, only: cross_product

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: occ_list(sys%nel)
        real(p), intent(out) :: psingle, pdouble

        integer :: i, j, virt_syms(2, sys%sym0_tot:sys%sym_max_tot), nsingles, ndoubles, isyma, isymb, ims1, ims2, isym1

        select case(sys%system)
        case(hub_k,ueg,ringium)
            ! Only double excitations
            psingle = 0.0_p
            pdouble = 1.0_p
        case(hub_real,heisenberg)
            ! Only single excitations
            psingle = 1.0_p
            pdouble = 0.0_p
        case(read_in)

            ! Count number of basis functions in each symmetry.
            virt_syms = sys%read_in%pg_sym%nbasis_sym_spin
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
                isym1 = sys%basis%basis_fns(occ_list(i))%sym
                nsingles = nsingles + virt_syms(ims1,isym1)
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
                        isymb = cross_product(sys, &
                                        sys%basis%basis_fns(occ_list(i))%sym, &
                                        sys%basis%basis_fns(occ_list(j))%sym)
                        if (sys%momentum_space) then
                            ! For momentum symmetry need to take care as symmetries are not self-inverse.
                            isymb = cross_product(sys, sys%read_in%mom_sym%inv_sym(isyma), isymb)
                        else
                            isymb = cross_product(sys, isyma, isymb)
                        end if

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

    subroutine find_parallel_spin_prob_mol(sys, pparallel)
        
        ! WARNING:
        ! Only call for read_in systems.

        ! Estimate pattempt_parallel by finding the ration of |Hij->ab|
        ! where ij are parallel to when they are not.

        ! In:
        !   sys: information about the system to be studied
        ! Out:
        !   pparallel: Estimate for pattempt_parallel

        use system, only: sys_t
        use proc_pointers, only: slater_condon2_excit_ptr, abs_hmatel_ptr
        use read_in_symmetry, only: cross_product_basis_read_in
        use hamiltonian_data, only: hmatel_t
#ifdef PARALLEL
        use parallel

        integer :: displs_nbasis(0:nprocs-1)
        integer :: sizes_nbasis(0:nprocs-1)
        integer :: ierr
        real(p) :: parallel_weight_tot, ortho_weight_tot
#endif

        type(sys_t), intent(in) :: sys
        real(p), intent(out) :: pparallel
        
        integer :: iproc_nbasis_start, iproc_nbasis_end
        type(hmatel_t)  :: hmatel
        integer :: i, j, a, b, i_tmp, j_tmp, a_tmp, b_tmp, ij_sym, isymb
        real(p) :: parallel_weight, ortho_weight

#ifdef PARALLEL
        ! Initialise do-loop range for each processor, [iproc_nbasis_start,iproc_nbasis_end].
        ! [todo] - get_proc_loop_range can also assign in serial mode. Disadvantage: would need to define displs_nbasis
        ! [todo] - and sizes_nbasis for serial mode.
        call get_proc_loop_range(sys%basis%nbasis, iproc_nbasis_start, iproc_nbasis_end, displs_nbasis, sizes_nbasis)

        parallel_weight_tot = 0.0_p
        ortho_weight_tot = 0.0_p
#else
        iproc_nbasis_start = 1
        iproc_nbasis_end = sys%basis%nbasis
#endif

        parallel_weight = 0.0_p
        ortho_weight = 0.0_p

        do i = iproc_nbasis_start, iproc_nbasis_end
            !$omp parallel do default(none) &
            !$omp shared(sys,i,slater_condon2_excit_ptr,abs_hmatel_ptr) &
            !$omp private(i_tmp, j_tmp, a_tmp, b_tmp, j, a, b, ij_sym, isymb, hmatel) reduction(+:parallel_weight,ortho_weight)
            do j = 1, sys%basis%nbasis
                if (i /= j) then
                    if (j < i) then
                        i_tmp = j
                        j_tmp = i
                    else
                        i_tmp = i
                        j_tmp = j
                    end if
                    ! The symmetry of b, isymb, is given by
                    ! (sym_i* x sym_j* x sym_a)* = sym_b
                    ! (at least for Abelian point groups)
                    ! ij_sym: symmetry conjugate of the irreducible representation spanned by the codensity
                    !        \phi_i*\phi_j. (We assume that ij is going to be in the bra of the excitation.)
                    ! [todo] - Check whether order of i and j matters here.
                    ij_sym = sys%read_in%sym_conj_ptr(sys%read_in, cross_product_basis_read_in(sys, i_tmp, j_tmp))
                    do a = 1, sys%basis%nbasis
                        if ((a /= i) .and. (a /= j)) then
                            isymb = sys%read_in%sym_conj_ptr(sys%read_in, &
                                        sys%read_in%cross_product_sym_ptr(sys%read_in, ij_sym, sys%basis%basis_fns(a)%sym))
                            do b = 1, sys%basis%nbasis
                                ! Check spin conservation and symmetry conservation.
                                if ((((sys%basis%basis_fns(i_tmp)%Ms == sys%basis%basis_fns(a)%Ms) .and. &
                                    (sys%basis%basis_fns(j_tmp)%Ms == sys%basis%basis_fns(b)%Ms)) .or. &
                                    ((sys%basis%basis_fns(i_tmp)%Ms == sys%basis%basis_fns(b)%Ms) .and. &
                                    (sys%basis%basis_fns(j_tmp)%Ms == sys%basis%basis_fns(a)%Ms))) .and. &
                                    (sys%basis%basis_fns(b)%sym == isymb) .and. ((b /= a) .and. (b /= i) .and. &
                                    (b /= j))) then
                                    if (b < a) then
                                        a_tmp = b
                                        b_tmp = a
                                    else
                                        a_tmp = a
                                        b_tmp = b
                                    end if
                                    hmatel = slater_condon2_excit_ptr(sys, i_tmp, j_tmp, a_tmp, b_tmp, .false.)
                                    if (sys%basis%basis_fns(i_tmp)%Ms == sys%basis%basis_fns(j_tmp)%Ms) then
                                        ! parallel spins
                                        parallel_weight = parallel_weight + abs_hmatel_ptr(hmatel)
                                    else
                                        ! not parallel
                                        ortho_weight = ortho_weight + abs_hmatel_ptr(hmatel)
                                    end if
                                end if
                            end do
                        end if
                    end do
                end if
            end do
            !$omp end parallel do
        end do

#ifdef PARALLEL
        call mpi_reduce(parallel_weight, parallel_weight_tot, 1, mpi_preal, MPI_SUM, root, MPI_COMM_WORLD, ierr)
        call mpi_reduce(ortho_weight, ortho_weight_tot, 1, mpi_preal, MPI_SUM, root, MPI_COMM_WORLD, ierr)
        if (iproc == root) pparallel = parallel_weight_tot/(parallel_weight_tot + ortho_weight_tot)
        call MPI_BCast(pparallel, 1, mpi_preal, root, MPI_COMM_WORLD, ierr)
#else
        pparallel = parallel_weight/(parallel_weight + ortho_weight)
#endif

    end subroutine find_parallel_spin_prob_mol

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
        real(dp) :: pextra
        integer :: nattempts

        nattempts = abs(int(population))
        pextra = abs(population) - nattempts
        ! If there is no probability of generating an extra attempt, then
        ! don't bother using an extra random number.
        if (abs(pextra) > depsilon) then
            if (pextra > get_rand_close_open(rng)) nattempts = nattempts + 1
        end if

    end function decide_nattempts

    subroutine load_balancing_report(nparticles, nstates_active, use_mpi_barriers, spawn_mpi_time, determ_mpi_time, io_unit)

        ! In:
        !    nparticles: number of particles in each space, on this process only.
        !    nstates_active: number of occupied states, on this process only.
        !    use mpi_barriers: if true then MPI_Barrier calls have been
        !        performed and timed, and their timings will be considered here.
        !    spawn_mpi_time: MPI timings for the spawned list, on this process
        !        only.
        ! In (optional):
        !    determ_mpi_time: MPI timings for semi-stochastic communications,
        !        on this process only.
        !    io_unit: io unit to write report to.

        ! Print out a load-balancing report when run in parallel showing how
        ! determinants and walkers/particles are distributed over the processors.

        use parallel
        use utils, only: int_fmt

        real(dp), intent(in) :: nparticles(:)
        integer, intent(in) :: nstates_active
        logical, intent(in) :: use_mpi_barriers
        type(parallel_timing_t), intent(in) :: spawn_mpi_time
        type(parallel_timing_t), optional, intent(in) :: determ_mpi_time
        integer, intent(in), optional :: io_unit
#ifdef PARALLEL
        real(dp) :: load_data(size(nparticles), nprocs)
        integer :: load_data_int(nprocs)
        integer :: i, ierr
        real(p) :: barrier_this_proc
        real(p) :: spawn_comms(nprocs), determ_comms(nprocs), barrier_time(nprocs)
        character(4) :: lfmt
        integer :: iunit

        iunit = 6
        if (present(io_unit)) iunit = io_unit

        if (nprocs > 1) then
            if (parent) then
                write (iunit,'(1X,a14,/,1X,14("^"),/)') 'Load balancing'
                write (iunit,'(1X,a77,/)') "The final distribution of walkers and determinants across the processors was:"
            endif
            call mpi_gather(nparticles, size(nparticles), mpi_real8, load_data, size(nparticles), &
                            mpi_real8, 0, MPI_COMM_WORLD, ierr)
            if (parent) then
                do i = 1, size(nparticles)
                    if (size(nparticles) > 1) write (iunit,'(1X,a,'//int_fmt(i,1)//')') 'Particle type:', i
                    write (iunit,'(1X,"Min # of particles on a processor:",6X,es13.6)') minval(load_data(i,:))
                    write (iunit,'(1X,"Max # of particles on a processor:",6X,es13.6)') maxval(load_data(i,:))
                    write (iunit,'(1X,"Mean # of particles on a processor:",5X,es13.6,/)') sum(load_data(i,:))/nprocs
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
                write (iunit,'(1X,"Min # of determinants on a processor:",3X,'//lfmt//')') minval(load_data_int)
                write (iunit,'(1X,"Max # of determinants on a processor:",3X,'//lfmt//')') maxval(load_data_int)
                write (iunit,'(1X,"Mean # of determinants on a processor:",2X,es13.6)') real(sum(load_data_int), p)/nprocs
                write (iunit,'()')
                if (use_mpi_barriers) then
                    write (iunit,'(1X,"Min time taken by MPI barrier calls:",5X,f8.2,"s")') minval(barrier_time)
                    write (iunit,'(1X,"Max time taken by MPI barrier calls:",5X,f8.2,"s")') maxval(barrier_time)
                    write (iunit,'(1X,"Mean time taken by MPI barrier calls:",4X,f8.2,"s")') sum(barrier_time)/nprocs
                    write (iunit,'()')
                end if
                write (iunit,'(1X,"Min time taken by walker communication:",5X,f8.2,"s")') minval(spawn_comms)
                write (iunit,'(1X,"Max time taken by walker communication:",5X,f8.2,"s")') maxval(spawn_comms)
                write (iunit,'(1X,"Mean time taken by walker communication:",4X,f8.2,"s")') sum(spawn_comms)/nprocs
                write (iunit,'()')
                if (present(determ_mpi_time)) then
                    write (iunit,'(1X,"Min time taken by semi-stochastic communication:",5X,f8.2,"s")') minval(determ_comms)
                    write (iunit,'(1X,"Max time taken by semi-stochastic communication:",5X,f8.2,"s")') maxval(determ_comms)
                    write (iunit,'(1X,"Mean time taken by semi-stochastic communication:",4X,f8.2,"s")') &
                        sum(determ_comms)/nprocs
                    write (iunit,'()')
                end if
            end if
        end if
#endif

    end subroutine load_balancing_report

    subroutine redistribute_particles(states, real_factor, pops, nstates, nparticles, spawn)

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
        ! In/Out:
        !    nparticles: number of excips on the current processor.
        !    pops: Population on occupied excitors.  On output the
        !        populations of excitors which are sent to other processors are
        !        set to zero.
        !    nstates: number of occupied excitors on the current processor.
        !    spawn: spawn_t object.  On output particles which need to be sent
        !        to another processor have been added to the correct position in
        !        the spawned store.

        ! Note: this adds particles to the spawn_t object which are not actually
        ! spawning events, so care must be taken with calculating the spawning
        ! rate if this procedure is used.

        use calc, only: dmqmc_calc, doing_calc
        use const, only: i0, dp
        use spawn_data, only: spawn_t
        use spawning, only: assign_particle_processor_dmqmc, assign_particle_processor, add_spawned_particles
        use parallel, only: iproc, nprocs

        integer(i0), intent(in) :: states(:,:)
        integer(int_p), intent(in) :: real_factor
        integer(int_p), intent(inout) :: pops(:,:)
        integer, intent(inout) :: nstates
        real(dp), intent(inout) :: nparticles(:)
        type(spawn_t), intent(inout) :: spawn

        real(p) :: nsent(size(nparticles))

        integer :: iexcitor, pproc, slot

        nsent = 0.0_dp

        !$omp parallel do default(none) &
        !$omp shared(nstates, states, pops, spawn, iproc, nprocs) &
        !$omp private(pproc, slot) reduction(+:nsent)
        do iexcitor = 1, nstates
            if (doing_calc(dmqmc_calc)) then
                ! NB this assumes sys%basis%info_string_len = 0 for DMQMC. This is
                ! for convenience of not changing any additional interfaces.
                ! If this is not the case in future this must be changed to be
                ! compatible.
                call assign_particle_processor_dmqmc(states(:,iexcitor), spawn%bit_str_nbits, 0, spawn%hash_seed, &
                                               spawn%hash_shift, spawn%move_freq, nprocs, pproc, slot, spawn%proc_map%map, &
                                               spawn%proc_map%nslots)
            else
                call assign_particle_processor(states(:,iexcitor), spawn%bit_str_nbits, spawn%hash_seed, spawn%hash_shift, &
                                               spawn%move_freq, nprocs, pproc, slot, spawn%proc_map%map, spawn%proc_map%nslots)
            end if
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

    subroutine redistribute_semi_stoch_t(sys, propagator, reference, annihilation_flags, psip_list, spawn, io_unit, determ)

        ! Recreate the semi_stoch_t object (if a non-empty space is in use).
        ! This requires sending deterministic states to their new processes
        ! and recreating the related objects, such as the deterministic
        ! Hamiltonian.

        ! In:
        !    sys: system being studied.
        !    propagator: propagator_t object containing quasinewton parameters
        !    reference: information on the current reference determinant.
        !    annihilation_flags: calculation specific annihilation flags.
        !    psip_list: list of particles and their locations.
        !    spawn: spawn_t object, required for determining the new processes
        !        labels for deterministic states.
        ! In/Out:
        !    determ: The deterministic space being used, as required for
        !        semi-stochastic calculations.

        use semi_stoch, only: dealloc_semi_stoch_t, init_semi_stoch_t

        use spawn_data, only: spawn_t
        use system, only: sys_t
        use qmc_data, only: semi_stoch_t, empty_determ_space, reuse_determ_space, particle_t
        use qmc_data, only: annihilation_flags_t, semi_stoch_in_t, propagator_t
        use reference_determinant, only: reference_t

        type(sys_t), intent(in) :: sys
        type(annihilation_flags_t), intent(in) :: annihilation_flags
        type(propagator_t), intent(in) :: propagator
        type(reference_t), intent(in) :: reference
        type(particle_t), intent(inout) :: psip_list
        type(spawn_t), intent(in) :: spawn
        type(semi_stoch_t), intent(inout) :: determ
        integer, intent(in) :: io_unit

        type(semi_stoch_in_t) :: ss_in_new

        if (determ%space_type /= empty_determ_space) then
            ! Create a temporary semi_stoch_in_t object to pass into the
            ! init_semi_stoch routine to update the determ object.
            ss_in_new%space_type = reuse_determ_space
            ss_in_new%projection_mode = determ%projection_mode

            ! Deallocate the semi_stoch_t instance, except for the list of all
            ! deterministic states, which we want to reuse.
            call dealloc_semi_stoch_t(determ, keep_dets=.true.)
            ! Recreate the semi_stoch_t instance, by reusing the deterministic
            ! space already generated, but with states on their new processes.
            call init_semi_stoch_t(determ, ss_in_new, sys, propagator, psip_list, reference, annihilation_flags, spawn, &
                                    .false., io_unit)
        end if

    end subroutine redistribute_semi_stoch_t

! --- Output routines ---

    subroutine initial_qmc_status(sys, qmc_in, qs, ntot_particles, doing_ccmc, io_unit)

        ! Calculate the projected energy based upon the initial walker
        ! distribution (either via a restart or as set during initialisation)
        ! and print out.

        ! In:
        !    sys: system being studied.
        !    qmc_in: input options relating to QMC methods.
        !    ntot_particles: total number of particles in each space.
        !    doing_ccmc: true if doing ccmc calculation.
        ! In/Out:
        !    qs: qmc_state_t object.
        ! In (optional):
        !    io_unit: io unit to write any reporting to.

        use parallel, only: parent
        use qmc_io, only: write_qmc_report
        use qmc_data, only: qmc_in_t, qmc_state_t
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        type(qmc_state_t), intent(inout), target :: qs
        real(dp), intent(in) :: ntot_particles(qs%psip_list%nspaces)
        logical, intent(in) :: doing_ccmc
        integer, optional, intent(in) :: io_unit
        
        if (parent) then
            if (doing_ccmc) then
                qs%estimators%nattempts = nint(qs%estimators%D0_population)
                call write_qmc_report(qmc_in, qs, 0, ntot_particles, 0.0, .true., .false., cmplx_est=sys%read_in%comp, &
                                        nattempts=.true., io_unit=io_unit)
            else
                call write_qmc_report(qmc_in, qs, 0, ntot_particles, 0.0, .true., .false., cmplx_est=sys%read_in%comp, &
                    io_unit=io_unit)
            end if
        end if

    end subroutine initial_qmc_status

    subroutine initial_ci_projected_energy(sys, qs, nb_comm, ntot_particles)

        ! Calculate the projected energy based upon the initial walker
        ! distribution.  proj_energy and D0_population are both accumulated in
        ! update_proj_energy.

        ! In:
        !    sys: system being studies
        !    qmc_in: input options relating to QMC methods.
        !    nb_comm: true if performing a calculation using non-blocking communications.
        ! In/Out:
        !    qs: qmc_state_t object. On output the estimator_t quantities are updated for each space based upon the particle
        !        distribution in the space.
        ! Out:
        !    ntot_particles: total number of particles in each space.

        use parallel

        use determinants, only: det_info_t, alloc_det_info_t, dealloc_det_info_t, decode_det
        use excitations, only: excit_t, get_excitation
        use importance_sampling, only: importance_sampling_weight
        use proc_pointers, only: update_proj_energy_ptr
        use qmc_data, only: qmc_state_t, nb_rep_t, zero_estimators_t
        use system, only: sys_t
        use hamiltonian_data

        type(sys_t), intent(in) :: sys
        type(qmc_state_t), intent(inout), target :: qs
        logical, intent(in) :: nb_comm
        real(dp), intent(out) :: ntot_particles(qs%psip_list%nspaces)

        integer :: idet, ispace
        real(p) :: real_population(qs%psip_list%nspaces), weighted_population(qs%psip_list%nspaces)
        type(det_info_t) :: cdet
        type(hmatel_t) :: hmatel
        type(excit_t) :: D0_excit
#ifdef PARALLEL
        integer :: ierr
        real(p) :: proj_energy_sum(qs%psip_list%nspaces), D0_population_sum(qs%psip_list%nspaces)
        complex(p) :: proj_energy_comp_sum(qs%psip_list%nspaces), D0_population_comp_sum(qs%psip_list%nspaces)
#endif
        call zero_estimators_t(qs%estimators)

        ! [todo] - HFS, DMQMC quantities
        call alloc_det_info_t(sys, cdet)
        do idet = 1, qs%psip_list%nstates
            cdet%f = qs%psip_list%states(:,idet)
            call decode_det(sys%basis, cdet%f, cdet%occ_list)
            cdet%data => qs%psip_list%dat(:,idet)
            real_population = real(qs%psip_list%pops(:,idet),p)/qs%psip_list%pop_real_factor
            do ispace = 1, qs%psip_list%nspaces
                weighted_population(ispace) = importance_sampling_weight(qs%trial, cdet, real_population(ispace))
            end do
            ! WARNING!  We assume only the bit string, occ list and data field
            ! are required to update the projected estimator.
            D0_excit = get_excitation(sys%nel, sys%basis, cdet%f, qs%ref%f0)
            if (sys%read_in%comp) then
                do ispace = 1, qs%psip_list%nspaces, 2
                    call update_proj_energy_ptr(sys, qs%ref%f0, qs%trial%wfn_dat, cdet, weighted_population(ispace:ispace+1), &
                                                qs%estimators(ispace), D0_excit, hmatel)
                    qs%estimators(ispace+1)%D0_population_comp = qs%estimators(ispace)%D0_population_comp
                    qs%estimators(ispace+1)%proj_energy_comp = qs%estimators(ispace)%proj_energy_comp
                end do
            else
                do ispace = 1, qs%psip_list%nspaces
                    call update_proj_energy_ptr(sys, qs%ref%f0, qs%trial%wfn_dat, cdet, [weighted_population(ispace)], &
                                                qs%estimators(ispace), D0_excit, hmatel)
                end do
            end if

        end do
        call dealloc_det_info_t(cdet)

#ifdef PARALLEL
        ! Non-blocking delays reporting by a report loop and initialisation of summation of the
        ! estimators is performed in init_non_blocking_comm.
        if (.not.nb_comm) then
            call mpi_allreduce(qs%estimators%proj_energy, proj_energy_sum, qs%psip_list%nspaces, mpi_preal, &
                               MPI_SUM, MPI_COMM_WORLD, ierr)
            call mpi_allreduce(qs%estimators%proj_energy_comp, proj_energy_comp_sum, qs%psip_list%nspaces, mpi_pcomplex, &
                               MPI_SUM, MPI_COMM_WORLD, ierr)
            call mpi_allreduce(qs%psip_list%nparticles, ntot_particles, qs%psip_list%nspaces, MPI_REAL8, &
                               MPI_SUM, MPI_COMM_WORLD, ierr)
            call mpi_allreduce(qs%estimators%D0_population, D0_population_sum, qs%psip_list%nspaces, mpi_preal, MPI_SUM, &
                               MPI_COMM_WORLD, ierr)
            call mpi_allreduce(qs%estimators%D0_population_comp, D0_population_comp_sum, qs%psip_list%nspaces, mpi_pcomplex, &
                               MPI_SUM, MPI_COMM_WORLD, ierr)
            call mpi_allreduce(qs%psip_list%nstates, qs%estimators%tot_nstates, qs%psip_list%nspaces, MPI_INTEGER, MPI_SUM, &
                               MPI_COMM_WORLD, ierr)
            qs%estimators%proj_energy = proj_energy_sum
            qs%estimators%proj_energy_comp = proj_energy_comp_sum
            qs%estimators%D0_population = D0_population_sum
            qs%estimators%D0_population_comp = D0_population_comp_sum
        end if
#else
        ntot_particles = qs%psip_list%nparticles
        qs%estimators%tot_nstates = qs%psip_list%nstates
#endif

    end subroutine initial_ci_projected_energy

    subroutine initial_cc_projected_energy(sys, qs, rng_seed, logging_info, cumulative_abs_real_pops, ntot_particles)

        ! Calculate the projected energy based upon the initial walker
        ! distribution for a CC wavefunction ansatz.

        ! In:
        !    sys: system being studies
        !    qmc_in: input options relating to QMC methods.
        !    rng_seed: seed for DSFMT random number generator. Use our own RNG stream so we don't disturb
        !        the subseuquent Markov chain in the QMC calculation.
        !    logging_info: logging status. Used only in debugging runs.
        ! In/Out:
        !    qs: qmc_state_t object. On output the estimator_t quantities are updated for each space based upon the particle
        !        distribution in the space.
        !    cumulative_abs_real_pops: scratch space for calculating the cumulative population distribution. Must be at least
        !        of size of current number of states. Do not use output unless the excip distribution has not been subseuquently
        !        modified!
        ! Out:
        !    ntot_particles: total number of particles in each space.

        use dSFMT_interface, only: dSFMT_t, dSFMT_init, dSFMT_end
        use parallel

        use ccmc_data, only: cluster_t, ex_lvl_dist_t
        use ccmc_utils, only: cumulative_population, get_D0_info
        use ccmc_selection, only: select_cluster, select_nc_cluster
        use determinants, only: det_info_t, alloc_det_info_t, dealloc_det_info_t
        use excitations, only: excit_t, get_excitation
        use hamiltonian_data, only: hmatel_t
        use logging, only: logging_t
        use proc_pointers, only: update_proj_energy_ptr
        use qmc_data, only: qmc_state_t, zero_estimators_t
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        type(qmc_state_t), intent(inout), target :: qs
        integer, intent(in) :: rng_seed
        type(logging_t), intent(in) :: logging_info
        real(p), intent(inout), allocatable :: cumulative_abs_real_pops(:)
        real(dp), intent(out) :: ntot_particles(qs%psip_list%nspaces)

        type(det_info_t) :: cdet
        type(cluster_t) :: cluster
        type(hmatel_t) :: hmatel
        type(excit_t) :: D0_excit
        real(p) :: tot_abs_real_pop, pop(2)
        type(ex_lvl_dist_t) :: ex_lvl_dist
        type(dSFMT_t) :: rng
        integer(int_64) :: iattempt, nattempts
        integer :: D0_pos, D0_proc, nD0_proc
        complex(p) :: D0_normalisation
#ifdef PARALLEL
        integer :: ierr
        real(p) :: proj_energy_sum(qs%psip_list%nspaces)
        complex(p) :: proj_energy_comp_sum(qs%psip_list%nspaces)
#endif

        call zero_estimators_t(qs%estimators)

        call dSFMT_init(rng_seed, 50000, rng)
        
        D0_pos = 1
        call get_D0_info(qs, sys%read_in%comp, D0_proc, D0_pos, nD0_proc, D0_normalisation)
        associate(pl=>qs%psip_list)
            call cumulative_population(pl%pops, pl%states(sys%basis%tot_string_len,:), pl%nstates, D0_proc, D0_pos, &
                                       pl%pop_real_factor, .false., sys%read_in%comp, cumulative_abs_real_pops, &
                                       tot_abs_real_pop, ex_lvl_dist)
        end associate
        ! Choose one cluster of size 2 for each excip we have.
        nattempts = nint(tot_abs_real_pop, int_64)

        ! Generate clusters of size 1 or 2; only these have a chance of connecting to the reference determinant.
        call alloc_det_info_t(sys, cdet)
        allocate(cluster%excitors(2))
        do iattempt = 1, nattempts + qs%psip_list%nstates
            if (iattempt < qs%psip_list%nstates) then
                call select_nc_cluster(sys, qs%psip_list, qs%ref%f0, iattempt, 0.0_p, .false., cdet, cluster)
            else
                ! Note: even if we're doing linked CC, the clusters contributing to the projected estimator must not contain
                ! excitors involving the same orbitals so we need only look for unlinked clusters.
                call select_cluster(rng, sys, qs%psip_list, qs%ref%f0, 2, .false., nattempts, D0_normalisation, 0.0_p, D0_pos, &
                                cumulative_abs_real_pops, tot_abs_real_pop, 2, 2, logging_info, cdet, cluster)
            end if
            if (cluster%excitation_level /= huge(0)) then

                D0_excit = get_excitation(sys%nel, sys%basis, cdet%f, qs%ref%f0)
                pop = [real(cluster%amplitude,p), aimag(cluster%amplitude)]*cluster%cluster_to_det_sign/cluster%pselect
                ! Note: replica tricks is not yet implemented in CCMC so only have (at most) real and imaginary spaces, which are
                ! handled in select_cluster/select_nc_cluster.
                call update_proj_energy_ptr(sys, qs%ref%f0, qs%trial%wfn_dat, cdet, pop, qs%estimators(1), D0_excit, hmatel)
            end if
        end do
        deallocate(cluster%excitors)
        call dealloc_det_info_t(cdet)
        call dSFMT_end(rng)
        ! WARNING: be careful when implementing ccmc replica or something else using nspaces!
        ! This is to be safe as get_D0_info sets D0_normalisation on all processors.
        qs%estimators(1)%D0_population = real(D0_normalisation, p)
        qs%estimators(1)%D0_population_comp = D0_normalisation
        if (sys%read_in%comp) then
            qs%estimators(2)%D0_population_comp = qs%estimators(1)%D0_population_comp
            qs%estimators(2)%proj_energy_comp = qs%estimators(1)%proj_energy_comp
        end if

#ifdef PARALLEL
        call mpi_allreduce(qs%estimators%proj_energy, proj_energy_sum, qs%psip_list%nspaces, mpi_preal, &
                           MPI_SUM, MPI_COMM_WORLD, ierr)
        call mpi_allreduce(qs%estimators%proj_energy_comp, proj_energy_comp_sum, qs%psip_list%nspaces, mpi_pcomplex, &
                           MPI_SUM, MPI_COMM_WORLD, ierr)
        call mpi_allreduce(qs%psip_list%nparticles, ntot_particles, qs%psip_list%nspaces, MPI_REAL8, &
                           MPI_SUM, MPI_COMM_WORLD, ierr)
        call mpi_allreduce(qs%psip_list%nstates, qs%estimators%tot_nstates, qs%psip_list%nspaces, MPI_INTEGER, MPI_SUM, &
                           MPI_COMM_WORLD, ierr)
        qs%estimators%proj_energy = proj_energy_sum
        qs%estimators%proj_energy_comp = proj_energy_comp_sum
#else
        ntot_particles = qs%psip_list%nparticles
        qs%estimators%tot_nstates = qs%psip_list%nstates
#endif
    
    end subroutine initial_cc_projected_energy

! --- QMC loop and cycle initialisation routines ---

    subroutine init_report_loop(qs, bloom_stats)

        ! Initialise a report loop (basically zero quantities accumulated over
        ! a report loop).

        use bloom_handler, only: bloom_stats_t, bloom_stats_init_report_loop
        use qmc_data, only: qmc_state_t, zero_estimators_t

        type(qmc_state_t), intent(inout) :: qs
        type(bloom_stats_t), intent(inout) :: bloom_stats

        call bloom_stats_init_report_loop(bloom_stats)

        ! Ensure D0_population from last cycle is set appropriately if restarting
        qs%estimators%D0_population_old = qs%estimators%D0_population

        qs%spawn_store%rspawn = 0.0_p
        call zero_estimators_t(qs%estimators)

    end subroutine init_report_loop

    subroutine init_mc_cycle(psip_list, spawn, nattempts, ndeath, min_attempts, complx)

        ! Initialise a Monte Carlo cycle (basically zero/reset cycle-level
        ! quantities).

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
        !    complx: true if using real and imaginary psips.

        use calc, only: doing_calc, ct_fciqmc_calc, ccmc_calc, dmqmc_calc
        use const, only: int_64, int_p
        use qmc_data, only: particle_t
        use spawn_data, only: spawn_t

        type(particle_t), intent(inout) :: psip_list
        type(spawn_t), intent(inout) :: spawn
        integer(int_64), intent(in), optional :: min_attempts
        logical, intent(in), optional :: complx
        integer(int_64), intent(out) :: nattempts
        integer(int_p), intent(out) :: ndeath

        logical :: complx_loc

        complx_loc = .false.
        if (present(complx)) complx_loc = complx

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
            ! (Note int is used rather than nint due to a minor error in implementing
            ! reals for CCMC so we are just maintaining behaviour.  The difference is
            ! really minimal...)
            nattempts = int(psip_list%nparticles(1), int_64)
            if (complx_loc) nattempts = nattempts + abs(nint(psip_list%nparticles(2), int_64))
        else if (doing_calc(dmqmc_calc)) then
            ! Each particle and each end gets to attempt to spawn onto a
            ! connected determinant and a chance to die/clone.
            nattempts = nint(4*psip_list%nparticles(1)*psip_list%nspaces, int_64)
        else
            ! Each particle gets to attempt to spawn onto a connected
            ! determinant and a chance to die/clone.
            nattempts = nint(2*psip_list%nparticles(1), int_64)
            if (complx_loc) nattempts = nattempts + abs(nint(2*psip_list%nparticles(2), int_64))
        end if

        if (present(min_attempts)) nattempts = max(nattempts, min_attempts)

    end subroutine init_mc_cycle

    subroutine load_balancing_wrapper(sys, propagator, reference, load_bal_in, annihilation_flags, nb_comm, io_unit, rng, &
                                      psip_list, spawn, par_info, determ)

        ! In:
        !    sys: system being studied
        !    propagator: propagotor_t containing quasinewton parameters
        !    reference: current reference determinant.
        !    load_bal_in: input options for load balancing.
        !    annihilation_flags: calculation specific annihilation flags.
        !    nb_comm: true if using non-blocking communications.
        ! In/Out:
        !    rng: random number generator.
        !    psip_list: total population (on this proccesor) is used to set
        !       nattempts and population is redistributed if requested by the
        !       load balancing approach.
        !    spawn: spawn_t object for holding spawned particles.  Reset on exit
        !       and updated with the new version of par_info%load%proc_map.
        !    par_info: type containing parallel information of the state of the
        !       system (load balancing and non-blocking).  Holds the 'master' copy
        !       of the updated proc_map on exit.
        ! In/Out (optional):
        !    determ: The deterministic space being used, as required for
        !        semi-stochastic calculations.

        ! WARNING: all spawn_t objects (bar the spawn argument provided) which operate
        ! on the same particle_t object must be updated with the proc_info from the
        ! master copy (par_info%load%proc_map).  It is the programmer's responsibility
        ! to ensure this happens.

        use system, only: sys_t
        use qmc_data, only: load_bal_in_t, annihilation_flags_t, particle_t
        use qmc_data, only: parallel_t, semi_stoch_t, propagator_t
        use spawn_data, only: spawn_t
        use dSFMT_interface, only: dSFMT_t
        use load_balancing, only: do_load_balancing
        use reference_determinant, only: reference_t

        type(sys_t), intent(in) :: sys
        type(propagator_t), intent(in) :: propagator
        type(reference_t), intent(in) :: reference
        type(load_bal_in_t), intent(in) :: load_bal_in
        type(annihilation_flags_t), intent(in) :: annihilation_flags
        logical, intent(in) :: nb_comm
        type(dSFMT_t), intent(inout) :: rng
        type(particle_t), intent(inout) :: psip_list
        type(spawn_t), intent(inout) :: spawn
        type(parallel_t), intent(inout) :: par_info
        type(semi_stoch_t), optional, intent(inout) :: determ
        integer, intent(in) :: io_unit

        if (par_info%load%needed) then
            call do_load_balancing(psip_list, spawn, par_info, load_bal_in)
            call redistribute_load_balancing_dets(rng, sys, propagator, reference, determ, psip_list, &
                                                  spawn, annihilation_flags, io_unit)
            ! If using non-blocking communications we still need this flag to
            ! be set.
            if (.not. nb_comm) par_info%load%needed = .false.
        end if

    end subroutine load_balancing_wrapper

! --- QMC loop and cycle termination routines ---

    subroutine end_report_loop(out_unit, qmc_in, iteration, update_tau, qs, ntot_particles, nspawn_events, semi_stoch_shift_it, &
                               semi_stoch_start_it, soft_exit, load_bal_in, update_estimators, bloom_stats, doing_lb, &
                               nb_comm, comp, error, vary_shift_reference)

        ! In:
        !    out_unit: File unit to write ouput to.
        !    qmc_in: input optons relating to QMC methods.
        !    iteration: The current iteration of the simulation.
        !    nspawn_events: The total number of spawning events to this process.
        !    semi_stoch_shift_it: How many iterations after the shift starts
        !        to vary to begin using semi-stochastic.
        !    load_bal_in: input options for load balancing.
        ! In/Out:
        !    update_tau: true if the processor thinks the timestep should be
        !        rescaled. This will be updated on output to only be true if
        !        in variable shift mode and if tau_search is being used.
        !    qs: qmc_state_t object. Energy estimators are updated.
        !    ntot_particles: total number (across all processors) of
        !        particles in the simulation at end of the previous report loop.
        !        Returns the current total number of particles for use in the
        !        next report loop if update_estimators is true.
        !    semi_stoch_start_it: The iteration on which to start performing
        !        semi-stochastic.
        ! Out:
        !    soft_exit: true if the user has requested an immediate exit of the
        !        QMC algorithm via the interactive functionality.
        ! In (optional):
        !    update_estimators: update the (FCIQMC/CCMC) energy estimators.  Default: true.
        !    doing_lb: true if doing load balancing.
        !    nb_comm: true if using non-blocking communications.
        !    load_bal_in: input options for load balancing.
        !    comp: true if doing qmc calculation with real and imaginary walkers.
        !    vary_shift_reference: if true, vary shift to control reference, not total, population
        ! In/Out (optional):
        !    bloom_stats: particle blooming statistics to accumulate.
        !    error: true if an error has occured and we need to quit.

        use energy_evaluation, only: update_energy_estimators, local_energy_estimators,         &
                                     update_energy_estimators_recv, update_energy_estimators_send
        use interact, only: calc_interact, check_interact, check_comms_file
        use parallel
        use system, only: sys_t
        use bloom_handler, only: bloom_stats_t, bloom_stats_warning
        use qmc_data, only: qmc_in_t, load_bal_in_t, qmc_state_t, nb_rep_t
        use spawning, only: update_pattempt

        integer, intent(in) :: out_unit
        type(qmc_in_t), intent(in) :: qmc_in
        integer, intent(in) :: iteration
        logical, intent(inout) :: update_tau
        type(qmc_state_t), intent(inout) :: qs
        integer, intent(in) :: nspawn_events
        logical, optional, intent(in) :: update_estimators, comp
        type(bloom_stats_t), optional, intent(inout) :: bloom_stats
        real(dp), intent(inout) :: ntot_particles(qs%psip_list%nspaces)
        integer, intent(in) :: semi_stoch_shift_it
        integer, intent(inout) :: semi_stoch_start_it
        logical, intent(out) :: soft_exit
        logical, intent(in), optional :: vary_shift_reference

        type(load_bal_in_t), intent(in) :: load_bal_in
        logical, optional, intent(in) :: doing_lb, nb_comm
        logical, optional, intent(inout) :: error

        logical :: update, vary_shift_before, nb_comm_local, comms_found, comp_param, overflow
        real(dp) :: rep_info_copy(size(qs%par_info%report_comm%rep_info))
        integer :: iunit

#ifdef PARALLEL
        integer :: ierr
#endif

        iunit = 6
        ! Only update the timestep if not in vary shift mode.
        update_tau = update_tau .and. .not. any(qs%vary_shift) .and. qmc_in%tau_search

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
        ! Check if complex parameter passed to function, and if not set to
        ! value that will have no effect.
        if (present(comp)) then
            comp_param = comp
        else
            comp_param = .false.
        end if
        if (present(update_estimators)) update = update_estimators
        if (update .and. .not. nb_comm_local) then
            call update_energy_estimators(qmc_in, qs, nspawn_events, ntot_particles, load_bal_in, doing_lb, &
                                          comms_found, error, update_tau, bloom_stats, vary_shift_reference, comp_param)
        else if (update) then
            ! Save current report loop quantitites.
            ! Can't overwrite the send buffer before message completion
            ! so copy information somewhere else.
            call local_energy_estimators(qs, rep_info_copy, nspawn_events, comms_found, error, update_tau, &
                                          bloom_stats, qs%par_info%report_comm%nb_spawn(2), comp_param)
            ! Receive previous iterations report loop quantities.
            call update_energy_estimators_recv(qmc_in, qs, qs%par_info%report_comm%request, ntot_particles, &
                                               qs%psip_list%nparticles_proc, load_bal_in, doing_lb, comms_found, error, &
                                               update_tau, bloom_stats, comp=comp_param)
            ! Send current report loop quantities.
            qs%par_info%report_comm%rep_info = rep_info_copy
            call update_energy_estimators_send(qs%par_info%report_comm)
        else
            update_tau = .false.
            call check_interact(comms_found)
        end if

        ! If we have just started varying the shift, then we can calculate the
        ! iteration at which to start semi-stochastic, if it is being turned on
        ! relative to the shift start.
        if ((.not. vary_shift_before) .and. all(qs%vary_shift) .and. (semi_stoch_shift_it /= -1)) &
            semi_stoch_start_it = semi_stoch_shift_it + iteration + 1

        if (qs%excit_gen_data%p_single_double%vary_psingles) then

#ifdef PARALLEL
            ! If any ps%rep_accum%overflow_loc is true (i.e. at least in one MPI proc there was a lack of precision in the
            ! number of single/double excitations), stop here and fix pattempt_single.
            call mpi_allreduce(qs%excit_gen_data%p_single_double%rep_accum%overflow_loc, overflow, 1, MPI_LOGICAL, MPI_LAND, &
                            MPI_COMM_WORLD, ierr)
#endif
            
            if ((qs%vary_shift(1)) .or. (overflow)) then
                if ((overflow) .and. (parent)) then
                    ! Make a note of the overflow.
                    write(iunit, '(1X, "# Had to stop varying pattempt_single due to lack of precision.")')
                end if
                ! Stop varying pattempt_single when the shift has started varying. This means that we do not update
                ! pattempt_single at the end of this report loop and of any future report loops.
                qs%excit_gen_data%p_single_double%vary_psingles = .false.
                ! Write (final) pattempt_single to output file.
                ! Format adapted from writing out shift damping.
                if (parent) write(iunit, '(1X, "# pattempt_single chosen to be:",1X,es17.10)') qs%excit_gen_data%pattempt_single
            else
                ! Possibly (if there were enough excitations) update pattempt_single.
                call update_pattempt(qs%excit_gen_data)
            end if
        end if
        
        call calc_interact(comms_found, out_unit, soft_exit, qs)

        if (qs%reblock_done) soft_exit = .true.

    end subroutine end_report_loop

    subroutine end_mc_cycle(nspawn_events, ndeath, real_factor, nattempts, rspawn)

        ! Execute common code at the end of a Monte Carlo cycle.

        ! In:
        !    nspawn_events: number of successful spawning events in the cycle.
        !    ndeath: (unscaled) number of particle deaths in the cycle.
        !    real_factor: the encoding factor by which the stored populations are multiplied
        !       to enable non-integer populations.
        !    nattempts: number of attempted spawning events in the cycle.
        ! In/Out:
        !    rspawn: running total of spawning rate.

        integer, intent(in) :: nspawn_events
        integer(int_p), intent(in) :: ndeath, real_factor
        integer(int_64), intent(in) :: nattempts
        real(p), intent(inout) :: rspawn

        ! Add the spawning rate (for the processor) to the running
        ! total.
        rspawn = rspawn + spawning_rate(nspawn_events, ndeath, real_factor, nattempts)

    end subroutine end_mc_cycle

    function spawning_rate(nspawn_events, ndeath, real_factor, nattempts) result(rate)

        ! Calculate the rate of spawning on the current processor.
        ! In:
        !    nspawn_events: number of successful spawning events during the
        !       MC cycle.
        !    ndeath: (unscaled) number of particles that were killed/cloned
        !       during the MC cycle.
        !    real_factor: The factor by which populations are multiplied to
        !        enable non-integer populations.
        !    nattempts: The number of attempts to spawn made in order to
        !       generate the current population of walkers in the spawned arrays.

        use const, only: p, int_p, int_64

        real(p) :: rate
        integer, intent(in) :: nspawn_events
        integer(int_p), intent(in) :: ndeath, real_factor
        integer(int_64), intent(in) :: nattempts
        real(p) :: ndeath_real

        ! Death is not scaled when using reals.
        ndeath_real = real(ndeath,p)/real_factor

        ! The total spawning rate is
        !   (nspawn + ndeath) / nattempts
        ! In, for example, the timestep algorithm each particle has 2 attempts
        ! (one to spawn on a different determinant and one to clone/die).
        ! ndeath is the number of particles that died, which hence equals the
        ! number of successful death attempts (assuming the timestep is not so
        ! large that death creates more than one particle).
        ! By ignoring the number of particles spawned in a single event, we
        ! hence are treating the death and spawning events on the same footing.
        if (nattempts > 0) then
            rate = (nspawn_events + ndeath_real)/nattempts
        else
            ! Can't have done anything.
            rate = 0.0_p
        end if

    end function spawning_rate

    subroutine rescale_tau(tau, factor)

        ! Scale the timestep by across all the processors.  This is performed if at least
        ! one processor thinks it should be.

        ! In:
        !    tau: timestep to be updated.
        !    factor: factor to scale the timestep by (default: 0.95).

        use parallel

        real(p), intent(inout) :: tau
        real(p), intent(in), optional :: factor
        integer :: iunit

        iunit = 6

        if (present(factor)) then
            tau = factor*tau
        else
            tau = 0.950_p*tau
        end if
        if (parent) write(iunit, '(1X, "# Warning timestep changed to:",1X,es17.10)') tau

    end subroutine rescale_tau

    subroutine redistribute_load_balancing_dets(rng, sys, propagator, reference, determ, &
                                                 psip_list, spawn, annihilation_flags, io_unit)

        ! When doing load balancing we need to redistribute chosen sections of
        ! main list to be sent to their new processors. This is a wrapper which
        ! takes care of this and resets the load balancing tag so that no
        ! further load balancing is attempted this report loop. Currently don't
        ! check if anything will actually be sent.

        ! Also if a non-empty semi-stochastic deterministic space is being used
        ! then this object needs to be redistributed, too.

        ! In:
        !    sys: system being studied.
        !    propagator: propagator_t giving quasinewton parameters
        !    reference: current reference determinant.
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
        !    spawn: spawn_t object.  Used to assign and send particles to their new
        !       processor.

        use annihilation, only: direct_annihilation
        use dSFMT_interface, only: dSFMT_t
        use qmc_data, only: semi_stoch_t, particle_t, annihilation_flags_t, propagator_t
        use reference_determinant, only: reference_t
        use spawn_data, only: spawn_t
        use system, only: sys_t

        type(dSFMT_t), intent(inout) :: rng
        type(sys_t), intent(in) :: sys
        type(propagator_t), intent(in) :: propagator
        type(reference_t), intent(in) :: reference
        type(semi_stoch_t), optional, intent(inout) :: determ
        type(particle_t), intent(inout) :: psip_list
        type(spawn_t), intent(inout) :: spawn
        type(annihilation_flags_t), intent(in) :: annihilation_flags
        integer, intent(in) :: io_unit

        associate(pl=>psip_list)
            call redistribute_particles(pl%states, pl%pop_real_factor, pl%pops, pl%nstates, pl%nparticles, spawn)
        end associate

        ! Merge determinants which have potentially moved processor back into
        ! the appropriate main list.
        call direct_annihilation(sys, rng, reference, annihilation_flags, psip_list, spawn)
        spawn%head = spawn%head_start

        if (present(determ)) call redistribute_semi_stoch_t(sys, propagator, reference, annihilation_flags, &
                                                             psip_list, spawn, io_unit, determ)

    end subroutine redistribute_load_balancing_dets

end module qmc_common
