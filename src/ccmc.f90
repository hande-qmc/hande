module ccmc

! Module  for performing coupled cluster Monte Carlo (CCMC) calculations).

! Due to the similarities with FCIQMC, we can re-use lots of the same routines
! (especially the spawning, death and annihilation).  As a result, the structure
! of do_ccmc is remarkably similar to the other do_*mc routines.

use const, only: p, lint

implicit none

contains

    subroutine do_ccmc()

        ! Run the CCMC algorithm starting from the initial walker distribution
        ! using the timestep algorithm.

        ! See notes about the implementation of this using function pointers
        ! in fciqmc_main.

        use parallel

        use annihilation, only: direct_annihilation
        use basis, only: basis_length, nbasis
        use calc, only: truncation_level
        use determinants, only: det_info, alloc_det_info, dealloc_det_info
        use excitations, only: excit, get_excitation_level
        use energy_evaluation, only: update_energy_estimators
        use fciqmc_data
        use fciqmc_common
        use fciqmc_restart, only: dump_restart
        use interact, only: fciqmc_interact
        use proc_pointers
        use spawning, only: create_spawned_particle_truncated

        integer :: ireport, icycle
        integer(lint) :: iattempt, nattempts, nparticles_old(sampling_size)
        type(det_info) :: cdet

        integer :: ndeath, nspawned, excitation_level
        real(p) :: cluster_amplitude, pcluster
        type(excit) :: connection

        logical :: soft_exit

        real :: t1, t2

        ! Allocate det_info components.
        call alloc_det_info(cdet)

        ! from restart
        nparticles_old = nparticles_old_restart

        ! Main fciqmc loop.
        if (parent) call write_fciqmc_report_header()
        call initial_fciqmc_status()
        ! Initialise timer.
        call cpu_time(t1)

        do ireport = 1, nreport

            ! Zero report cycle quantities.
            proj_energy = 0.0_p
            rspawn = 0.0_p
            D0_population = 0.0_p

            do icycle = 1, ncycles

                ! Reset the current position in the spawning array to be the
                ! slot preceding the first slot.
                spawning_head = spawning_block_start

                ! Number of spawning attempts that will be made.
                ! Each particle gets to attempt to spawn onto a connected
                ! determinant and a chance to die/clone.
                ! This is used for accounting later, not for controlling the spawning.
                nattempts = 2*nparticles(1)

                ! Reset death counter
                ndeath = 0

                ! Allow one spawning & death attempt for each excip on the
                ! processor.
                do iattempt = 1, nparticles(1)

                    ! TODO: initiator.
                    call select_cluster(nparticles(1), cdet, excitation_level, cluster_amplitude, pcluster)

                    if (excitation_level <= truncation_level+2) then

                        ! TODO: projected estimator.
                        call update_proj_energy_ptr(cdet, cluster_amplitude)

                        call spawner_ccmc(cdet, cluster_amplitude, pcluster, nspawned, connection)

                        if (nspawned /= 0) then
                            ! TODO: initiator (use create_spawned_particle_ptr
                            ! proc pointer).
                            call create_spawned_particle_truncated(cdet, connection, nspawned, spawned_pop)
                        end if

                        ! Does the cluster collapsed onto D0 produce
                        ! a determinant is in the truncation space?  If so, also
                        ! need to attempt a death/cloning step.
                        if (excitation_level <= truncation_level) call stochastic_ccmc_death(cdet, cluster_amplitude, pcluster)

                    end if

                end do

                ! Add the spawning rate (for the processor) to the running
                ! total.
                rspawn = rspawn + spawning_rate(ndeath, nattempts)

                ! D0_population is communicated in the direct_annihilation
                ! algorithm for efficiency.
                call direct_annihilation()

            end do

            ! Update the energy estimators (shift & projected energy).
            call update_energy_estimators(nparticles_old)

            call cpu_time(t2)

            ! t1 was the time at the previous iteration, t2 the current time.
            ! t2-t1 is thus the time taken by this report loop.
            if (parent) call write_fciqmc_report(ireport, nparticles_old(1), t2-t1)

            ! cpu_time outputs an elapsed time, so update the reference timer.
            t1 = t2

            call fciqmc_interact(soft_exit)
            if (soft_exit) exit
            if (mod(ireport, select_ref_det_every_nreports) == 0) call select_ref_det()

        end do

        if (parent) then
            call write_fciqmc_final(ireport)
            write (6,'()')
        end if

        call load_balancing_report()

        if (soft_exit) then
            mc_cycles_done = mc_cycles_done + ncycles*ireport
        else
            mc_cycles_done = mc_cycles_done + ncycles*nreport
        end if

        if (dump_restart_file) call dump_restart(mc_cycles_done, nparticles_old(1))

        call dealloc_det_info(cdet)

    end subroutine do_ccmc

    subroutine select_cluster(nattempts, cdet, excitation_level, amplitude, pcluster)

        ! TODO: interface comments

        use calc, only: truncation_level
        use determinants, only: det_info
        use dSFMT_interface, only: genrand_real2
        use fciqmc_data, only: f0, D0_population, tot_walkers, &
                               walker_population
        use proc_pointers, only: decoder_ptr
        use utils, only: factorial

        integer(lint), intent(in) :: nattempts
        type(det_info), intent(inout) :: cdet
        integer, intent(out) :: excitation_level
        real(p), intent(out) :: amplitude, pcluster

        real(p) :: rand, psize
        integer :: cluster_size, i
        logical :: valid_cluster
        integer :: excitor(truncation_level+1)

        ! We shall accumulate the factors which comprise pcluster as we go.
        !   pcluster = n_sel p_size p_clust
        ! where
        !   n_sel   is the number of cluster selections made (by this
        !           processor?);
        !   p_size  is the probability of choosing a cluster of that size;
        !   p_clust is the probability of choosing a specific cluster given
        !           the choice of size.

        pcluster = nattempts

        ! Select the cluster size, i.e. the number of excitors in a cluster.
        ! For a given truncation_level, only clusters containing at most
        ! truncation_level+1 excitors.
        ! Following the process described by Thom in 'Initiator Stochastic
        ! Coupled Cluster Theory' (unpublished), each size, n_s, has probability
        ! p(n_s) = 1/2^(n_s+1), n_s=0,truncation_level and p(truncation_level+1)
        ! is such that \sum_{n_s=0}^{truncation_level+1} p(n_s) = 1.
        rand = genrand_real2()
        psize = 0.0_p
        cluster_size = -1
        do i = 0, truncation_level
            psize = psize + 1.0_p/2**(i+1)
            if (rand < psize) then
                ! Found size!
                cluster_size = i
                pcluster = pcluster/2**(i+1)
                exit
            end if
        end do
        ! If not set, then must be the largest possible cluster
        if (cluster_size < 0) then
            cluster_size = truncation_level + 1
            pcluster = pcluster/(1.0_p - psize)
        end if

        valid_cluster = .true.
        select case(cluster_size)
        case(0)
            ! Must be the reference.
            cdet%f = f0
            excitation_level = 0
            amplitude = D0_population ! TODO: don't accumulate D0 over report loops
            ! Only one cluster of this size to choose => p_clust = 1
        case default
            ! Select cluster from the excitors on the current processor in
            ! a uniform manner.
            amplitude = 1.0_p
            do i = 1, cluster_size
                ! Select a position between in the excitors list.
                excitor(i) = int(genrand_real2()*tot_walkers) + 1
                amplitude = amplitude*walker_population(1,excitor(i))
            end do
            ! Normalisation factor for amplitudes...
            amplitude = amplitude/(D0_population**(cluster_size-1)) ! TODO: don't accumulate D0 over report loops

            ! TODO: collapse cluster (return also excitation_level)

            ! We chose excitors uniformly.
            !  -> prob. of each excitor is 1./(number of excitors on processor).
            ! Overall probability is the product of this times the number of
            ! ways the cluster could have been selected, ie
            !   n_s!/n_excitors^n_s
            pcluster = (pcluster*factorial(cluster_size))/(tot_walkers**cluster_size)
        end select

        ! Fill in information about the cluster if required.
        if (valid_cluster) then
            call decoder_ptr(cdet%f, cdet)
        else
            ! Simply set excitation level to a too high (fake) level to avoid
            ! this cluster being used.
            excitation_level = huge(0)
        end if

    end subroutine select_cluster

    subroutine spawner_ccmc(cdet, amplitude, pcluster, nspawn, connection)

        ! Attempt to spawn a new particle on a connected excitor with
        ! probability
        !     \tau |<D'|H|D_s> A_s|
        !   -------------------------
        !   n_sel p_s p_clust p_excit
        ! where |D_s> is the determinant formed by applying the excitor to the
        ! reference determinant, A_s is the amplitude and D' is the determinant
        ! formed from applying a connected excitor to the reference determinant.
        ! See comments in select_cluster about n_sel, p_s and p_clust.  p_excit
        ! is the probability of choosing D' given D_s.

        ! This is just a thin wrapper around a system-specific excitation
        ! generator and a utility function.  We need to modify the spawning
        ! probability compared to the FCIQMC algorithm as we spawn from multiple
        ! excips at once (in FCIQMC we allow each psip to spawn individually)
        ! and have additional probabilities to take into account.

        ! In:
        !    cdet: info on the current excitor (cdet) that we will spawn
        !        from.
        !    amplitude: amplitude of cluster.
        !    pcluster: Overall probabilites of selecting this cluster, ie
        !        n_sel.p_s.p_clust.
        ! Out:
        !    nspawn: number of particles spawned.  0 indicates the spawning
        !        attempt was unsuccessful.
        !    connection: excitation connection between the current excitor
        !        and the child excitor, on which progeny are spawned.

        use determinants, only: det_info
        use excitations, only: excit
        use proc_pointers, only: gen_excit_ptr
        use spawning, only: attempt_to_spawn

        type(det_info), intent(in) :: cdet
        real(p), intent(in) :: amplitude
        real(p), intent(in) :: pcluster
        integer, intent(out) :: nspawn
        type(excit), intent(out) :: connection

        ! We incorporate the sign of the amplitude into the Hamiltonian matrix
        ! element, so we 'pretend' to attempt_to_spawn that all excips are
        ! actually spawned by positive excips.
        integer, parameter :: parent_sign = 1
        real(p) :: hmatel, pgen

        ! 1. Generate random excitation.
        call gen_excit_ptr(cdet, pgen, connection, hmatel)

        ! 2, Apply additional factors.
        hmatel = hmatel*amplitude
        pgen = pgen/pcluster

        ! 3. Attempt spawning.
        nspawn = attempt_to_spawn(hmatel, pgen, parent_sign)

    end subroutine spawner_ccmc

    subroutine stochastic_ccmc_death(cdet, amplitude, pcluster)

        ! Attempt to 'die' (ie create an excip on the current exictor, cdet%f)
        ! with probability
        !    \tau |<D_s|H|D_s> A_s|
        !    ----------------------
        !       n_sel p_s p_clust
        ! where |D_s> is the determinant formed by applying the excitor to the
        ! reference determinant and A_s is the amplitude.  See comments in
        ! select_cluster about the probabilities.

        ! In:
        !    cdet: info on the current excitor (cdet) that we will spawn
        !        from.
        !    amplitude: amplitude of cluster.
        !    pcluster: Overall probabilites of selecting this cluster, ie
        !        n_sel.p_s.p_clust.

        use determinants, only: det_info
        use fciqmc_data, only: tau, shift, spawned_pop, H00
        use excitations, only: excit
        use proc_pointers, only: sc0_ptr, create_spawned_particle_ptr
        use dSFMT_interface, only: genrand_real2

        type(det_info), intent(in) :: cdet
        real(p), intent(in) :: amplitude
        real(p), intent(in) :: pcluster

        real(p) :: pdeath, Kii
        integer :: nkill
        type(excit), parameter :: null_excit = excit( 0, [0,0,0,0], [0,0,0,0], .false.)

        ! TODO: optimise for the case where the cluster is either the reference
        ! determinant or consisting of a single excitor.
        Kii = sc0_ptr(cdet%f) - H00 - shift

        pdeath = (tau*Kii*amplitude)/pcluster

        ! Number that will definitely die
        ! Create nkill excips with opposite sign to parent excip...
        nkill = -int(pdeath)

        ! Stochastic death...
        pdeath = pdeath - nkill
        if (abs(pdeath) > genrand_real2()) then
            ! Increase magnitude of nkill...
            nkill = nkill + sign(1, nkill)
        end if

        ! The excitor might be a composite cluster so we'll just create
        ! excips in the spawned list and allow the annihilation process to take
        ! care of the rest.
        ! Pass through a null excitation so that we create a spawned particle on
        ! the current excitor.
        call create_spawned_particle_ptr(cdet, null_excit, nkill, spawned_pop)

    end subroutine stochastic_ccmc_death

end module ccmc
