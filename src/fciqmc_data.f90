module fciqmc_data

! Data for fciqmc calculations and procedures which manipulate fciqmc and only
! fciqmc data.

use const
use csr, only: csrp_t
use dmqmc_data, only: rdm_t
use spawn_data, only: spawn_t
use hash_table, only: hash_table_t
use calc, only: parallel_t
use parallel, only: parallel_timing_t
implicit none

!--- Input data: FCIQMC ---

! [todo] - It is somewhat inelegant to store/pass around real_factor or real_bit_shift separately,
! [todo] - when really they are variables telling us about the representation of the population data.
! [todo] - As such, the separation is not entrely helpful.  It would make more sense encoded with the
! [todo] - data itself -  it makes no sense without the data, and the data make no sense without it.
! [todo] - Convert to a derived type.  (From AJWT.)
!
! [todo] - The code currently uses c = real(a)/real_factor everywhere.  However, one of 
! [todo] - the joys of fixed precision arithmetic is that c=(a*b)>>real_bit_shift.
! [todo] - is (or used to be) a lot faster than c=(real(a)/real_factor)*(real(b)/real_factor)
! [todo] - It is however a little more obfuscated though, but carrying fixed-precision all the
! [todo] - through the code might give a benefit.  Of course the best implementation would
! [todo] - be to hide the details in a class.  (From AJWT.)
!
! [todo] - real_bit_shift and real_factor really should be compile-time constants.
! Real amplitudes can be any multiple of 2**(-real_bit_shift). They are
! encoded as integers by multiplying them by 2**(real_bit_shift).
integer :: real_bit_shift
! real_factor = 2**(real_bit_shift)
integer(int_p) :: real_factor

!--- Walker data ---

! Rate of spawning.  This is a running total over MC cycles on each processor
! until it is summed over processors and averaged over cycles in
! update_energy_estimators.
real(p) :: rspawn

! In DMQMC the trace of the density matrix is an important quantity
! used in calculating all thermal estimators. This quantity stores
! the this value, Tr(\rho), where rho is the density matrix which
! the DMQMC algorithm calculates stochastically.
real(p), allocatable :: trace(:) ! (walker_global%nspaces)

! The following indicies are used to access components of DMQMC numerators.
enum, bind(c)
    enumerator :: energy_ind = 1
    enumerator :: energy_squared_ind
    enumerator :: correlation_fn_ind
    enumerator :: staggered_mag_ind
    enumerator :: full_r2_ind
    enumerator :: terminator ! unused except in num_dmqmc_operators
   ! NOTE: if you add a new estimator then you must insert it before terminator.
end enum

! This variable holds the total number of operators which are implemented
! for DMQMC.
integer, parameter :: num_dmqmc_operators = terminator - 1

! numerators stores the numerators for the estimators in DMQMC. These
! are, for a general operator O which we wish to find the thermal average of:
! \sum_{i,j} \rho_{ij} * O_{ji}
! This variabe will store this value from the first iteration of each
! report loop. At the end of a report loop, the values from each
! processor are combined and stored in numerators on the parent
! processor. This is then output, and the values of numerators
! are reset on each processor to start the next report loop.
real(p) :: numerators(num_dmqmc_operators)

! When using the replica_tricks option, if the rdm in the first
! simulation if denoted \rho^1 and the ancillary rdm is denoted
! \rho^2 then renyi_2 holds:
! x = \sum_{ij} \rho^1_{ij} * \rho^2_{ij}.
! The indices of renyi_2 hold this value for the various rdms being
! calculated. After post-processing averaging, this quantity should
! be normalised by the product of the corresponding RDM traces.
! call it y. Then the renyi-2 entropy is then given by -log_2(x/y).
real(p), allocatable :: renyi_2(:)

! rdm_traces(i,j) holds the trace of replica i of the rdm with label j.
real(p), allocatable :: rdm_traces(:,:) ! (walker_global%nspaces, nrdms)

! If this logical is true then the program runs the DMQMC algorithm with
! importance sampling.
! dmqmc_sampling_prob stores the factors by which the probabilities of
! spawning to a larger excitation are reduced by. So, when spawning from
! a diagonal element to a element with one excitation, the probability
! of spawning is reduced by a factor sampling_probs(1).
! accumulated_probs(i) stores the multiplication of all the elements
! of sampling_probs up to the ith element. This quantity is often
! needed, so it is stored.
real(p), allocatable :: accumulated_probs(:) ! (min(nel, nsites-nel) + 1)
real(p), allocatable :: accumulated_probs_old(:) ! (min(nel, nsites-nel) + 1)

real(dp), allocatable :: weight_altering_factors(:)

! Calculate replicas (ie evolve two wavefunctions/density matrices at once)?
! Currently only implemented for DMQMC.
logical :: replica_tricks = .false.

real(p), allocatable :: excit_dist(:) ! (0:max_number_excitations)

! This stores the reduces matrix, which is slowly accumulated over time
! (on each processor).
real(p), allocatable :: reduced_density_matrix(:,:)

! Spawned lists for rdms.
type rdm_spawn_t
    type(spawn_t) :: spawn
    ! Spawn with the help of a hash table to avoid a sort (which is extremely
    ! expensive when a large number of keys are repeated--seem to hit worst case
    ! performance in quicksort).
    type(hash_table_t) :: ht
end type rdm_spawn_t
type(rdm_spawn_t), allocatable :: rdm_spawn(:)

! The total number of rdms beings calculated.
integer :: nrdms = 0

! This stores  information for the various RDMs that the user asks to be
! calculated. Each element of this array corresponds to one of these RDMs.
type(rdm_t), allocatable :: fci_rdm_info(:) ! (nrdms)

! The total number of translational symmetry vectors.
! This is only set and used when performing rdm calculations.
integer :: nsym_vec

! The unit of the file reduced_dm.
integer :: rdm_unit

! When using the Neel singlet trial wavefunction, it is convenient
! to store all possible amplitudes in the wavefunction, since
! there are relativley few of them and they are expensive to calculate
real(dp), allocatable :: neel_singlet_amp(:) ! (nsites/2) + 1

!--- Restart data ---

! Restart data.
integer :: mc_cycles_done = 0

! Type for storing parallel information: see calc for description.
type(parallel_t) :: par_info

contains

    !--- Statistics. ---

    function spawning_rate(nspawn_events, ndeath, nattempts) result(rate)

        ! Calculate the rate of spawning on the current processor.
        ! In:
        !    nspawn_events: number of successful spawning events during the
        !       MC cycle.
        !    ndeath: (unscaled) number of particles that were killed/cloned
        !       during the MC cycle.
        !    nattempts: The number of attempts to spawn made in order to
        !       generate the current population of walkers in the spawned arrays.

        use parallel, only: nprocs

        real(p) :: rate
        integer, intent(in) :: nspawn_events
        integer(int_p), intent(in) :: ndeath
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

    !--- Output procedures ---

    subroutine write_fciqmc_report_header(ntypes, dmqmc_in)

        ! In:
        !    ntypes: number of particle types being sampled.
        !    dmqmc_in (optional): input options relating to DMQMC.

        use calc, only: doing_calc, hfs_fciqmc_calc, dmqmc_calc, doing_dmqmc_calc
        use calc, only: dmqmc_energy, dmqmc_energy_squared, dmqmc_staggered_magnetisation
        use calc, only: dmqmc_correlation, dmqmc_full_r2, dmqmc_rdm_r2
        use dmqmc_data, only: dmqmc_in_t, dmqmc_estimates_global
        use utils, only: int_fmt

        integer, intent(in) :: ntypes
        type(dmqmc_in_t), optional, intent(in) :: dmqmc_in

        integer :: i, j
        character(16) :: excit_header

        write (6,'()')

        if (doing_calc(dmqmc_calc)) then
           write (6,'(1X,a12,3X,a13,17X,a5)', advance = 'no') &
           '# iterations','Instant shift','Trace'

            if (doing_dmqmc_calc(dmqmc_full_r2)) then
                write (6, '(13X,a7,14X,a7)', advance = 'no') 'Trace 2','Full S2'
            end if
            if (doing_dmqmc_calc(dmqmc_energy)) then
                write (6, '(2X,a19)', advance = 'no') '\sum\rho_{ij}H_{ji}'
            end if
            if (doing_dmqmc_calc(dmqmc_energy_squared)) then
                write (6, '(2X,a19)', advance = 'no') '\sum\rho_{ij}H2{ji}'
            end if
            if (doing_dmqmc_calc(dmqmc_correlation)) then
                write (6, '(2X,a19)', advance = 'no') '\sum\rho_{ij}S_{ji}'
            end if
            if (doing_dmqmc_calc(dmqmc_staggered_magnetisation)) then
                write (6, '(2X,a19)', advance = 'no') '\sum\rho_{ij}M2{ji}'
            end if
            if (doing_dmqmc_calc(dmqmc_rdm_r2)) then
                do i = 1, nrdms
                    write (6, '(16X,a3,'//int_fmt(i,0)//',1x,a2)', advance = 'no') 'RDM', i, 'S2'
                end do
            end if
            if (dmqmc_in%rdm%calc_inst_rdm) then
                do i = 1, nrdms
                    do j = 1, ntypes
                        write (6, '(7X,a3,'//int_fmt(i,0)//',1x,a5,1x,'//int_fmt(j,0)//')', advance = 'no') &
                                'RDM', i, 'trace', j
                    end do
                end do
            end if
            if (present(dmqmc_in)) then
                if (dmqmc_in%calc_excit_dist) then
                    do i = 0, ubound(dmqmc_estimates_global%excit_dist,1)
                        write (excit_header, '("Excit. level",1X,'//int_fmt(i,0)//')') i
                        write (6, '(5X,a16)', advance='no') excit_header
                    end do
                end if
            end if

            write (6, '(3X,a11,6X)', advance='no') '# particles'

        else
            write (6,'(1X,a13,3(2X,a17))', advance='no') &
                     "# iterations ", "Shift            ", "\sum H_0j N_j    ", "N_0              "
            if (doing_calc(hfs_fciqmc_calc)) then
                write (6,'(6(2X,a17))', advance='no') &
                    "H.F. Shift       ","\sum O_0j N_j    ","\sum H_0j N'_j   ","N'_0             ", &
                    "# H psips        ","# HF psips       "
            else
                write (6,'(4X,a9,8X)', advance='no') "# H psips"
            end if
        end if
        write (6,'(3X,"# states  # spawn_events  R_spawn   time")')

    end subroutine write_fciqmc_report_header

    subroutine write_fciqmc_report(qmc_in, qs, ireport, ntot_particles, elapsed_time, comment, non_blocking_comm, dmqmc_in)

        ! Write the report line at the end of a report loop.

        ! In:
        !    qmc_in: input options relating to QMC methods.
        !    qs: QMC state (containing shift and various estimators).
        !    ireport: index of the report loop.
        !    ntot_particles: total number of particles in main walker list.
        !    elapsed_time: time taken for the report loop.
        !    comment: if true, then prefix the line with a #.
        !    non_blocking_comm: true if using non-blocking communications
        !    dmqmc_in: input options relating to DMQMC.

        use calc, only: doing_calc, dmqmc_calc, hfs_fciqmc_calc, doing_dmqmc_calc
        use calc, only: dmqmc_energy, dmqmc_energy_squared, dmqmc_full_r2, dmqmc_rdm_r2
        use calc, only: dmqmc_correlation, dmqmc_staggered_magnetisation
        use dmqmc_data, only: dmqmc_in_t
        use qmc_data, only: qmc_in_t, walker_global, qmc_state_t
        use dmqmc_data, only: dmqmc_estimates_global

        type(qmc_in_t), intent(in) :: qmc_in
        type(qmc_state_t), intent(in) :: qs
        integer, intent(in) :: ireport
        real(p), intent(in) :: ntot_particles(:)
        real, intent(in) :: elapsed_time
        logical, intent(in) :: comment, non_blocking_comm
        type(dmqmc_in_t), optional, intent(in) :: dmqmc_in

        integer :: mc_cycles, i, j, ntypes

        ntypes = size(ntot_particles)

        ! For non-blocking communications we print out the nth report loop
        ! after the (n+1)st iteration. Adjust mc_cycles accordingly
        if (.not. non_blocking_comm) then
            mc_cycles = ireport*qmc_in%ncycles
        else
            mc_cycles = (ireport-1)*qmc_in%ncycles
        end if

        if (comment) then
            write (6,'(1X,"#",1X)', advance='no')
        else
            write (6,'(3X)', advance='no')
        end if

        ! See also the format used in inital_fciqmc_status if this is changed.

        ! DMQMC output.
        if (doing_calc(dmqmc_calc)) then
            write (6,'(i10,2X,es17.10,2X,es17.10)',advance = 'no') &
                (mc_cycles_done+mc_cycles-qmc_in%ncycles), qs%shift(1), dmqmc_estimates_global%trace(1)
            ! The trace on the second replica.
            if (doing_dmqmc_calc(dmqmc_full_r2)) then
                write(6, '(3X,es17.10)',advance = 'no') dmqmc_estimates_global%trace(2)
            end if

            ! Renyi-2 entropy for the full density matrix.
            if (doing_dmqmc_calc(dmqmc_full_r2)) then
                write (6, '(4X,es17.10)', advance = 'no') dmqmc_estimates_global%numerators(full_r2_ind)
            end if

            ! Energy.
            if (doing_dmqmc_calc(dmqmc_energy)) then
                write (6, '(4X,es17.10)', advance = 'no') dmqmc_estimates_global%numerators(energy_ind)
            end if

            ! Energy squared.
            if (doing_dmqmc_calc(dmqmc_energy_squared)) then
                write (6, '(4X,es17.10)', advance = 'no') dmqmc_estimates_global%numerators(energy_squared_ind)
            end if

            ! Correlation function.
            if (doing_dmqmc_calc(dmqmc_correlation)) then
                write (6, '(4X,es17.10)', advance = 'no') dmqmc_estimates_global%numerators(correlation_fn_ind)
            end if

            ! Staggered magnetisation.
            if (doing_dmqmc_calc(dmqmc_staggered_magnetisation)) then
                write (6, '(4X,es17.10)', advance = 'no') dmqmc_estimates_global%numerators(staggered_mag_ind)
            end if

            ! Renyi-2 entropy for all RDMs being sampled.
            if (doing_dmqmc_calc(dmqmc_rdm_r2)) then
                do i = 1, nrdms
                    write (6, '(6X,es17.10)', advance = 'no') renyi_2(i)
                end do
            end if

            ! Traces for instantaneous RDM estimates.
            if (dmqmc_in%rdm%calc_inst_rdm) then
                do i = 1, nrdms
                    do j = 1, ntypes
                        write (6, '(2x,es17.10)', advance = 'no') rdm_traces(j,i)
                    end do
                end do
            end if

            ! The distribution of walkers on different excitation levels of the
            ! density matrix.
            if (present(dmqmc_in)) then
                if (dmqmc_in%calc_excit_dist) then
                    dmqmc_estimates_global%excit_dist = dmqmc_estimates_global%excit_dist/ntot_particles(1)
                    do i = 0, ubound(dmqmc_estimates_global%excit_dist,1)
                        write (6, '(4X,es17.10)', advance = 'no') dmqmc_estimates_global%excit_dist(i)
                    end do
                end if
            end if

            write (6, '(2X,es17.10)', advance='no') ntot_particles(1)

        else if (doing_calc(hfs_fciqmc_calc)) then
            write (6,'(i10,2X,6(es17.10,2X),es17.10,4X,es17.10,X,es17.10)', advance = 'no') &
                                             mc_cycles_done+mc_cycles, qs%shift(1),   &
                                             qs%estimators%proj_energy, qs%estimators%D0_population, &
                                             qs%shift(2), qs%estimators%proj_hf_O_hpsip, qs%estimators%proj_hf_H_hfpsip, &
                                             qs%estimators%D0_hf_population, &
                                             ntot_particles
        else
            write (6,'(i10,2X,2(es17.10,2X),es17.10,4X,es17.10)', advance='no') &
                                             mc_cycles_done+mc_cycles, qs%shift(1),   &
                                             qs%estimators%proj_energy, qs%estimators%D0_population, &
                                             ntot_particles
        end if
        write (6,'(2X,i10,4X,i12,2X,f7.4,2X,f6.3)') qs%estimators%tot_nstates, qs%estimators%tot_nspawn_events, rspawn, elapsed_time/qmc_in%ncycles

    end subroutine write_fciqmc_report

    subroutine end_fciqmc(nb_comm, reference, psip_list, spawn)

        ! Deallocate fciqmc data arrays.

        ! In:
        !    nb_comm: true if using non-blocking communications.
        ! In/Out (optional):
        !   reference: reference state. On exit, allocatable components are deallocated.
        !   psip_list: main particle_t object.  On exit, allocatable components are deallocated.
        !   spawn: spawn_t object.  On exit, allocatable components are deallocated.

        use checking, only: check_deallocate
        use spawn_data, only: dealloc_spawn_t
        use calc, only: dealloc_parallel_t
        use qmc_data, only: reference_t, particle_t
        use reference_determinant, only: dealloc_reference_t

        logical, intent(in) :: nb_comm
        type(reference_t), intent(inout), optional :: reference
        type(particle_t), intent(inout), optional :: psip_list
        type(spawn_t), intent(inout), optional :: spawn

        integer :: ierr

        if (present(reference)) then
            call dealloc_reference_t(reference)
        end if
        if (present(psip_list)) then
            if (allocated(psip_list%nparticles)) then
                deallocate(psip_list%nparticles, stat=ierr)
                call check_deallocate('psip_list%nparticles',ierr)
            end if
            if (allocated(psip_list%states)) then
                deallocate(psip_list%states, stat=ierr)
                call check_deallocate('psip_list%states',ierr)
            end if
            if (allocated(psip_list%pops)) then
                deallocate(psip_list%pops, stat=ierr)
                call check_deallocate('psip_list%pops',ierr)
            end if
            if (allocated(psip_list%dat)) then
                deallocate(psip_list%dat, stat=ierr)
                call check_deallocate('psip_list%dat',ierr)
            end if
            if (allocated(neel_singlet_amp)) then
                deallocate(neel_singlet_amp, stat=ierr)
                call check_deallocate('neel_singlet_amp',ierr)
            end if
            if (allocated(psip_list%nparticles_proc)) then
                deallocate(psip_list%nparticles_proc, stat=ierr)
                call check_deallocate('psip_list%nparticles_proc', ierr)
            end if
        end if
        call dealloc_parallel_t(nb_comm, par_info)
        if (present(spawn)) call dealloc_spawn_t(spawn)

    end subroutine end_fciqmc

end module fciqmc_data
