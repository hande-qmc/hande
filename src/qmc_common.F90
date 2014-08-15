module qmc_common

! Module containing routines common to different qmc algorithms.

use fciqmc_data

implicit none

public
private :: stochastic_round_int_32, stochastic_round_int_64

interface stochastic_round
    module procedure stochastic_round_int_32
    module procedure stochastic_round_int_64
end interface stochastic_round

contains

! --- Utility routines ---

    subroutine select_ref_det(sys)

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

        use calc, only: doing_calc, hfs_fciqmc_calc
        use determinants, only: decode_det, write_det
        use system, only: sys_t

        use parallel
        use errors, only: stop_all

        type(sys_t), intent(in) :: sys

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

        allocate(fmax(lbound(walker_dets, dim=1):ubound(walker_dets, dim=1)))

        H00_old = H00

        updated = .false.
        ! Find determinant with largest population.
        max_pop = 0_int_p
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

        if (all(fmax == f0)) then
            ! Max population on this processor is already the reference.  Don't change.
            in_data = (/ 0_int_p, int(iproc,int_p) /)
        else if (abs(max_pop) > ref_det_factor*abs(D0_population)) then
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
            f0 = fmax
            H00 = H00_max
            ! Broadcast updated data
            call mpi_bcast(f0, size(f0), mpi_det_integer, D0_proc, MPI_COMM_WORLD, ierr)
            call mpi_bcast(H00, 1, mpi_preal, D0_proc, MPI_COMM_WORLD, ierr)
        end if

#else

        if (abs(max_pop) > ref_det_factor*abs(D0_population) .and. any(fmax /= f0)) then
            updated = .true.
            f0 = fmax
            H00 = H00_max
        end if

#endif

        if (updated) then
            ! Update occ_list.
            call decode_det(sys%basis, f0, occ_list0)
            ! walker_data(1,i) holds <D_i|H|D_i> - H00_old.  Update.
            ! H00 is currently <D_0|H|D_0> - H00_old.
            ! Want walker_data(1,i) to be <D_i|H|D_i> - <D_0|H|D_0>
            ! We'll fix H00 later and avoid an extra tot_walkers*additions.
            do i = 1, tot_walkers
                walker_data(1,i) = walker_data(1,i) - H00
            end do
            ! The fold line in the folded spectrum approach is set (during
            ! initialisation) relative to the reference.
            fold_line = fold_line - H00
            ! Now set H00 = <D_0|H|D_0> so that future references to it are
            ! correct.
            H00 = H00 + H00_old
            if (doing_calc(hfs_fciqmc_calc)) call stop_all('select_ref_det', 'Not implemented for HFS.')
            if (parent) then
                write (6,'(1X,"#",1X,62("-"))')
                write (6,'(1X,"#",1X,"Changed reference det to:",1X)',advance='no')
                call write_det(sys%basis, sys%nel, f0, new_line=.true.)
                write (6,'(1X,"#",1X,"Population on old reference det (averaged over report loop):",f10.2)') D0_population
                write (6,'(1X,"#",1X,"Population on new reference det:",27X,i8)') max_pop
                write (6,'(1X,"#",1X,"E0 = <D0|H|D0> = ",f20.12)') H00
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

    subroutine cumulative_population(pops, nactive, D0_proc, D0_pos, cumulative_pops, tot_pop)

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
        !    D0_pos: position in the pops list of the reference.  Only relevant if
        !       1<=D0_pos<=nactive and the processor holds the reference.
        ! Out:
        !    cumulative_pops: running total of excitor population, i.e.
        !        cumulative_pops(i) = sum(abs(pops(1:i))), excluding the
        !        population on the reference if appropriate.
        !    tot_pop: total population (possibly excluding the population on the
        !       reference).

        ! NOTE: currently only the populations in the first psip/excip space are
        ! considered.  This should be changed if we do multiple simulations at
        ! once/Hellmann-Feynman sampling/etc.

        ! WARNING: almost certainly not suitable for a parallel implementation.

        use parallel, only: iproc

        integer(int_p), intent(in) :: pops(:,:)
        integer, intent(in) :: nactive, D0_proc, D0_pos
        integer(int_p), intent(out) :: cumulative_pops(:), tot_pop

        integer :: i

        cumulative_pops(1) = abs(pops(1,1))
        if (D0_proc == iproc) then
            ! Let's be a bit faster: unroll loops and skip over the reference
            ! between the loops.
            do i = 2, d0_pos-1
                cumulative_pops(i) = cumulative_pops(i-1) + abs(pops(1,i))
            end do
            ! Set cumulative on the reference to be the running total merely so we
            ! can continue accessing the running total from the i-1 element in the
            ! loop over excitors in slots above the reference.
            if (d0_pos == 1) cumulative_pops(d0_pos) = 0
            if (d0_pos > 1) cumulative_pops(d0_pos) = cumulative_pops(d0_pos-1)
            do i = d0_pos+1, nactive
                cumulative_pops(i) = cumulative_pops(i-1) + abs(pops(1,i))
            end do
        else
            ! V simple on other processors: no reference to get in the way!
            do i = 2, nactive
                cumulative_pops(i) = cumulative_pops(i-1) + abs(pops(1,i))
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
        real(dp), intent(in) :: population
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

    subroutine stochastic_round_int_32(rng, population, cutoff, ntypes)

        ! For any values in population less than cutoff, round up to cutoff or
        ! down to zero. This is done such that the expectation value of the
        ! resulting populations is equal to the input values.

        ! In/Out:
        !    rng: random number generator.
        ! In:
        !    population: populations to be stochastically rounded.
        !    cutoff: the value to round up to.
        !    ntypes: the number of values in population to apply this op to.

        use const, only: int_4
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open

        type(dSFMT_t), intent(inout) :: rng
        integer(int_4), intent(inout) :: population(:)
        integer(int_4), intent(in) :: cutoff
        integer, intent(in) :: ntypes
        integer :: itype
        real(p) :: r

        do itype = 1, ntypes
            if (abs(population(itype)) < cutoff .and. population(itype) /= 0_int_4) then
                r = get_rand_close_open(rng)*cutoff
                if (abs(population(itype)) > r) then
                    population(itype) = sign(cutoff, population(itype))
                else
                    population(itype) = 0_int_4
                end if
            end if
        end do

    end subroutine stochastic_round_int_32

    subroutine stochastic_round_int_64(rng, population, cutoff, ntypes)

        ! For any values in population less than cutoff, round up to cutoff or
        ! down to zero. This is done such that the expectation value of the
        ! resulting populations is equal to the input values.

        ! In/Out:
        !    rng: random number generator.
        ! In:
        !    population: populations to be stochastically rounded.
        !    cutoff: the value to round up to.
        !    ntypes: the number of values in population to apply this op to.

        use const, only: int_8
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open

        type(dSFMT_t), intent(inout) :: rng
        integer(int_8), intent(inout) :: population(:)
        integer(int_8), intent(in) :: cutoff
        integer, intent(in) :: ntypes
        integer :: itype
        real(p) :: r

        do itype = 1, ntypes
            if (abs(population(itype)) < cutoff .and. population(itype) /= 0_int_8) then
                r = get_rand_close_open(rng)*cutoff
                if (abs(population(itype)) > r) then
                    population(itype) = sign(cutoff, population(itype))
                else
                    population(itype) = 0_int_8
                end if
            end if
        end do

    end subroutine stochastic_round_int_64

    subroutine load_balancing_report()

        ! Print out a load-balancing report when run in parallel showing how
        ! determinants and walkers/particles are distributed over the processors.

#ifdef PARALLEL
        use parallel
        use utils, only: int_fmt

        real(dp) :: load_data(sampling_size, nprocs)
        integer(lint) :: load_data_lint(nprocs)
        integer :: i, ierr
        real(dp) :: comms(nprocs)
        character(4) :: lfmt

        if (nprocs > 1) then
            if (parent) then
                write (6,'(1X,a14,/,1X,14("^"),/)') 'Load balancing'
                write (6,'(1X,a77,/)') "The final distribution of walkers and determinants across the processors was:"
            endif
            call mpi_gather(nparticles, sampling_size, mpi_real8, load_data, sampling_size, &
                            mpi_real8, 0, MPI_COMM_WORLD, ierr)
            if (parent) then
                do i = 1, sampling_size
                    if (sampling_size > 1) write (6,'(1X,a,'//int_fmt(i,1)//')') 'Particle type:', i
                    write (6,'(1X,"Min # of particles on a processor:",6X,es12.6)') minval(load_data(i,:))
                    write (6,'(1X,"Max # of particles on a processor:",6X,es12.6)') maxval(load_data(i,:))
                    write (6,'(1X,"Mean # of particles on a processor:",5X,es12.6,/)') real(sum(load_data(i,:)), p)/nprocs
                end do
            end if
            call mpi_gather(tot_walkers, 1, mpi_integer8, load_data_lint, 1, mpi_integer8, 0, MPI_COMM_WORLD, ierr)
            call mpi_gather(annihilation_comms_time, 1, mpi_real8, comms, 1, mpi_real8, 0, MPI_COMM_WORLD, ierr)
            if (parent) then
                lfmt = int_fmt(maxval(load_data_lint),0)
                write (6,'(1X,"Min # of determinants on a processor:",3X,'//lfmt//')') minval(load_data_lint)
                write (6,'(1X,"Max # of determinants on a processor:",3X,'//lfmt//')') maxval(load_data_lint)
                write (6,'(1X,"Mean # of determinants on a processor:",2X,es12.6)') real(sum(load_data_lint), p)/nprocs
                write (6,'()')
                write (6,'(1X,"Min time taken by walker communication:",5X,f8.2,"s.")') minval(comms)
                write (6,'(1X,"Max time taken by walker communication:",5X,f8.2,"s.")') maxval(comms)
                write (6,'(1X,"Mean time taken by walker communication:",4X,f8.2,"s.")') real(sum(comms), p)/nprocs
                write (6,'()')
            end if
        end if
#endif

    end subroutine load_balancing_report

! --- Output routines ---

    subroutine initial_fciqmc_status(sys)

        ! Calculate the projected energy based upon the initial walker
        ! distribution (either via a restart or as set during initialisation)
        ! and print out.

        ! In:
        !    sys: system being studied.

        use determinants, only: det_info_t, alloc_det_info_t, dealloc_det_info_t, decode_det
        use excitations, only: excit
        use parallel
        use proc_pointers, only: update_proj_energy_ptr
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer :: idet
        real(dp) :: ntot_particles(sampling_size)
        real(dp) :: real_population(sampling_size)
        type(det_info_t) :: cdet
        real(p) :: hmatel
        type(excit) :: D0_excit
#ifdef PARALLEL
        integer :: ierr
        real(p) :: proj_energy_sum
#endif

        ! Calculate the projected energy based upon the initial walker
        ! distribution.  proj_energy and D0_population_cycle are both accumulated in
        ! update_proj_energy.
        proj_energy = 0.0_p
        call alloc_det_info_t(sys, cdet)
        do idet = 1, tot_walkers
            cdet%f = walker_dets(:,idet)
            call decode_det(sys%basis, cdet%f, cdet%occ_list)
            cdet%data => walker_data(:,idet)
            real_population = real(walker_population(:,idet),dp)/real_factor
            ! WARNING!  We assume only the bit string, occ list and data field
            ! are required to update the projected estimator.
            call update_proj_energy_ptr(sys, f0, cdet, real_population(1), &
                                        D0_population_cycle, proj_energy, D0_excit, hmatel)
        end do
        call dealloc_det_info_t(cdet)

#ifdef PARALLEL
        call mpi_allreduce(proj_energy, proj_energy_sum, sampling_size, mpi_preal, MPI_SUM, MPI_COMM_WORLD, ierr)
        proj_energy = proj_energy_sum
        call mpi_allreduce(nparticles, ntot_particles, sampling_size, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        call mpi_allreduce(D0_population_cycle, D0_population, 1, mpi_preal, MPI_SUM, MPI_COMM_WORLD, ierr)
        ! TODO: HFS, DMQMC quantities
#else
        ntot_particles = nparticles
        D0_population = D0_population_cycle
#endif

        if (parent) then
            ! See also the format used in write_fciqmc_report if this is changed.
            ! We prepend a # to make it easy to skip this point when do data
            ! analysis.
            call write_fciqmc_report(0, ntot_particles, 0.0, .true.)
        end if

    end subroutine initial_fciqmc_status

! --- QMC loop and cycle initialisation routines ---

    subroutine init_report_loop(bloom_stats)

        ! Initialise a report loop (basically zero quantities accumulated over
        ! a report loop).

        use bloom_handler, only: bloom_stats_t

        type(bloom_stats_t), intent(inout) :: bloom_stats

        bloom_stats%nwarnings_curr = 0

        proj_energy = 0.0_p
        rspawn = 0.0_p
        D0_population = 0.0_p

        ! DMQMC-specific...
        if (calculate_excit_distribution) excit_distribution = 0.0_p
        if (allocated(trace)) trace = 0.0_p
        if (allocated(estimator_numerators)) estimator_numerators = 0.0_p

    end subroutine init_report_loop

    subroutine init_mc_cycle(nattempts, ndeath, min_attempts)

        ! Initialise a Monte Carlo cycle (basically zero/reset cycle-level
        ! quantities).

        ! In:
        !    min_attempts (optional): if present, set nattempts to be at least this value.
        ! Out:
        !    nattempts: number of spawning attempts to be made (on the current
        !        processor) this cycle.
        !    ndeath: number of particle deaths that occur in a Monte Carlo
        !        cycle.  Reset to 0 on output.

        use calc, only: doing_calc, ct_fciqmc_calc, ccmc_calc, dmqmc_calc

        integer(lint), intent(in), optional :: min_attempts
        integer(lint), intent(out) :: nattempts
        integer(int_p), intent(out) :: ndeath

        ! Reset the current position in the spawning array to be the
        ! slot preceding the first slot.
        qmc_spawn%head = qmc_spawn%head_start

        ! Reset death counter
        ndeath = 0_int_p

        ! Reset accumulation of the reference population over the MC cycle.
        ! Reset accumulation of the projected estimator over MC cycle.
        ! This is only relevant in CCMC, where the estimators must be weighted
        ! by the number of times they're sampled.  (In contrast, these are
        ! currently calculated exactly in FCIQMC.)
        proj_energy_cycle = 0.0_p
        D0_population_cycle = 0.0_p

        ! Number of spawning attempts that will be made.
        ! For FCIQMC, this is used for accounting later, not for controlling the
        ! spawning.
        if (doing_calc(ct_fciqmc_calc) .or. doing_calc(ccmc_calc)) then
            ! ct algorithm: kinda poorly defined.
            ! ccmc: number of excitor clusters we'll randomly generate and
            ! attempt to spawn from.
            nattempts = nparticles(1)
        else if (doing_calc(dmqmc_calc)) then
            ! Each particle and each end gets to attempt to spawn onto a
            ! connected determinant and a chance to die/clone.
            nattempts = nint(4*nparticles(1)*sampling_size)
        else
            ! Each particle gets to attempt to spawn onto a connected
            ! determinant and a chance to die/clone.
            nattempts = nint(2*nparticles(1))
        end if

        if (present(min_attempts)) nattempts = max(nattempts, min_attempts)

    end subroutine init_mc_cycle

! --- QMC loop and cycle termination routines ---

    subroutine end_report_loop(sys, ireport, update_tau, ntot_particles, report_time, soft_exit, update_estimators)

        ! In:
        !    sys: system being studied.
        !    ireport: index of current report loop.
        !    update_tau: true if the processor thinks the timestep should be rescaled.
        !             Only used if not in variable shift mode and if tau_search is being
        !             used.
        !    update_estimators (optional): update the (FCIQMC/CCMC) energy estimators.  Default: true.
        ! In/Out:
        !    ntot_particles: total number (across all processors) of
        !        particles in the simulation at end of the previous report loop.
        !        Returns the current total number of particles for use in the
        !        next report loop if update_estimators is true.
        !    report_time: time at the start of the current report loop.  Returns
        !        the current time (ie the time for the start of the next report
        !        loop.
        ! Out:
        !    soft_exit: true if the user has requested an immediate exit of the
        !        QMC algorithm via the interactive functionality.

        use energy_evaluation, only: update_energy_estimators
        use interact, only: fciqmc_interact
        use parallel, only: parent
        use restart_hdf5, only: dump_restart_hdf5, restart_info_global, restart_info_global_shift
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: ireport
        logical, intent(in) :: update_tau
        logical, optional, intent(in) :: update_estimators
        real(dp), intent(inout) :: ntot_particles(sampling_size)
        real, intent(inout) :: report_time
        logical, intent(out) :: soft_exit

        real :: curr_time
        logical :: update

        ! Update the energy estimators (shift & projected energy).
        update = .true.
        if (present(update_estimators)) update = update_estimators
        if (update) call update_energy_estimators(ntot_particles)

        call cpu_time(curr_time)

        ! report_time was the time at the previous iteration.
        ! curr_time - report_time is thus the time taken by this report loop.
        if (parent) call write_fciqmc_report(ireport, ntot_particles, curr_time-report_time, .false.)

        ! Write restart file if required.
        if (dump_restart_file_shift .and. vary_shift) then
            dump_restart_file_shift = .false.
            call dump_restart_hdf5(restart_info_global_shift, mc_cycles_done+ncycles*ireport, ntot_particles)
        else if (mod(ireport,restart_info_global%write_restart_freq) == 0) then
            call dump_restart_hdf5(restart_info_global, mc_cycles_done+ncycles*ireport, ntot_particles)
        end if
        ! cpu_time outputs an elapsed time, so update the reference timer.
        report_time = curr_time

        call fciqmc_interact(soft_exit)
        if (.not.soft_exit .and. mod(ireport, select_ref_det_every_nreports) == 0) call select_ref_det(sys)

        if (tau_search .and. .not. vary_shift) call send_tau(update_tau)

    end subroutine end_report_loop

    subroutine end_mc_cycle(nspawn_events, ndeath, nattempts)

        ! Execute common code at the end of a Monte Carlo cycle.

        ! In:
        !    nspawn_events: number of successful spawning events in the cycle.
        !    ndeath: (unscaled) number of particle deaths in the cycle.
        !    nattempts: number of attempted spawning events in the cycle.

        integer, intent(in) :: nspawn_events
        integer(int_p), intent(in) :: ndeath
        integer(lint), intent(in) :: nattempts

        ! Add the spawning rate (for the processor) to the running
        ! total.
        rspawn = rspawn + spawning_rate(nspawn_events, ndeath, nattempts)

        ! Accumulate population on reference and projected estimator.
        D0_population = D0_population + D0_population_cycle
        proj_energy = proj_energy + proj_energy_cycle

    end subroutine end_mc_cycle

    subroutine send_tau(update_tau, factor)

        ! Scale the timestep by across all the processors.  This is performed if at least
        ! one processor thinks it should be.

        ! In:
        !    update_tau: true if the current processor thinks the timestep should be scaled.
        !    factor: factor to scale the timestep by (default: 0.95).

        use parallel

        integer :: ierr
        logical, intent(in) :: update_tau
        real(p), intent(in), optional :: factor

        logical :: update_tau_global

#ifdef PARALLEL
       call mpi_allreduce(update_tau, update_tau_global, 1, mpi_logical, MPI_LOR, MPI_COMM_WORLD, ierr)
#else
        update_tau_global = update_tau
#endif
        if(update_tau_global) then
            if (present(factor)) then
                tau = factor*tau
            else
                tau = 0.950_p*tau
            end if
            if(parent) write(6, '(1X, "# Warning timestep changed to: ",f7.5)'), tau
        end if

    end subroutine send_tau

end module qmc_common
