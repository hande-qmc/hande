module qmc_data

use const
use spawn_data, only: spawn_t, proc_map_t
use csr, only: csrp_t
use parallel, only: parallel_timing_t

implicit none

! --- QMC input ---

type qmc_in_t

    ! Seed used to initialise the dSFMT random number generator.
    integer :: seed

    ! True if allowing non-integer values for psip populations.
    logical :: real_amplitudes = .false.
    ! The minimum amplitude of a spawning event which can be added to
    ! the spawned list.
    ! If real amplitudes are not used then the following default will be
    ! overwritten by 0.0_p. In this case it will effectively not be used and all
    ! spawnings events will be integers.
    real(p) :: spawn_cutoff = 0.01_p

    ! Don't bother renormalising generation probabilities; instead allow forbidden
    ! excitations to be generated and then rejected.
    logical :: no_renorm = .false.

    ! probability of attempting single or double excitations...
    ! set to be nonsense value so can easily detect if it's given as an input option
    real(p) :: pattempt_single = -1, pattempt_double = -1

    ! timestep
    ! Note: qmc_state_t%tau is used and (if desired) updated during the course of a simulation)
    real(p) :: tau
    ! Are we doing a timestep search
    logical :: tau_search = .false.

    ! The shift (energy offset) is set to vary_shift_from when variable shift mode is entered.
    ! WARNING: if both initial_shift and vary_shift_from are set, then we expect the
    ! user to have been sensible.
    real(p) :: vary_shift_from = 0.0_p
    ! If true, then the when variable shift mode is entered the shift is set to be the
    ! current projected energy estimator.  Overrides vary_shift_from.
    logical :: vary_shift_from_proje = .false.
    ! Initial shift.
    real(p) :: initial_shift = 0.0_p
    ! Factor by which the changes in the population are damped when updating the
    ! shift.
    real(p) :: shift_damping = 0.050_dp

    ! Array sizes: main and spawned particle lists.
    ! If these are < 0, then the values represent the number of MB to be used.
    ! CARE: as we don't modify qmc_in_t objects, one should inspect the sizes
    ! used in particle_t and spawned_particle_t for the exact values used (which may
    ! be rounded for various reasons).
    integer :: walker_length
    integer :: spawned_walker_length

    ! The initial population on the reference determinant/trace of the density matrix.
    ! Overridden by a restart file.
    real(p) :: D0_population
    ! Number of particles before which varyshift mode is turned on.
    real(p) :: target_particles = huge(1.0_p)

    ! Using the initiator approximation?
    logical :: initiator_approx = .false.
    ! Population above which a determinant is an initiator.
    real(p) :: initiator_pop = 3.0_p

    ! number of monte carlo cycles/report loop
    integer :: ncycles
    ! number of report cycles
    ! the shift is updated and the calculation information printed out
    ! at the end of each report cycle.
    integer :: nreport

end type qmc_in_t

type fciqmc_in_t

    ! How often do we change the reference determinant to the determinant with
    ! greatest population?
    ! Default: we don't.
    integer :: select_ref_det_every_nreports = huge(1)

    ! Also start with D0_population on i_s|D_0>, where i_s is the spin-version
    ! operator.  This is only done if no restart file is used *and* |D_0> is not
    ! a closed shell determinant.
    logical :: init_spin_inv_D0 = .false.

    ! Factor by which the population on a determinant must exceed the reference
    ! determinant's population in order to be accepted as the new reference
    ! determinant.
    real(p) :: ref_det_factor = 1.50_p

    ! Flag for using non-blocking communications.
    ! Default: False.
    logical :: non_blocking_comm = .false.
    ! Flag for using load balancing.
    ! Default: False.
    logical :: doing_load_balancing = .false.

    ! Importance sampling (see enumerators below for allowed values).  Currently
    ! only relevant to the Heisenberg model.  trial_function will always be
    ! single_basis for other models to represent a single determinant.
    integer :: trial_function = 0
    ! If we are not using importance sampling, this is set to no_guiding, otherwise
    ! to a specific enumerator below to specify the corresponding guiding function
    ! being used.
    integer :: guiding_function = 0

end type fciqmc_in_t

type semi_stoch_in_t
    ! The iteration on which to turn on the semi-stochastic algorithm using the
    ! deterministic space specified by space_type.
    integer :: start_iter = 1
    ! The iteration on which to turn on the semi-stochastic algorithm, relative
    ! to the iteration that the shift starts to vary (or the iteration at which
    ! all shifts have started to vary, in the case of multiple replicas).
    integer :: shift_iter = -1
    ! space_type is used to tell the semi-stochastic initialisation routine
    ! which type of deterministic space to use. See the 'determ-space' parameters
    ! defined in semi_stoch.F90 for the various values it can take.
    integer :: space_type = 0
    ! Certain deterministic space types need a target size to be input to tell the
    ! semi-stochastic initialisation routine how many states to try and include. In
    ! such cases this variable should be set on input.
    integer :: target_size = 0
    ! If true then the deterministic states will be written to a file.
    logical :: write_determ_space = .false.
    ! If true then deterministic spawnings will not be added to the spawning list
    ! but rather treated separately via an extra MPI call.
    logical :: separate_annihil = .true.
end type semi_stoch_in_t

type ccmc_in_t
    ! How frequently (in log_2) an excitor can be moved to a different processor.
    ! See comments in spawn_t and assign_particle_processor.
    integer :: move_freq = 5
    ! Value of cluster%amplitude/cluster%pselect above which spawns are split up
    ! The default value corresponds to off.
    real(p) :: cluster_multispawn_threshold = huge(1.0_p)
    ! Use the full non-composite algorithm?
    logical :: full_nc = .false.
    ! Sample only linked clusters in CCMC?
    logical :: linked = .false.
end type ccmc_in_t

type restart_in_t
    ! Restart calculation from file.
    logical :: read_restart = .false.
    ! Index to read from (huge indicates not set).
    integer :: read_id = huge(0)
    ! Print out restart file.
    logical :: write_restart = .false.
    ! Index to write to (huge indicates not set).
    integer :: write_id = huge(0)
    ! Print out restart file every X iterations (in addition to at the end).
    integer :: write_freq = huge(0)
    ! Print out a restart file just before the shift turns on.
    logical :: write_restart_shift = .false.
    ! Index to write to (huge indicates not set).
    integer :: write_shift_id = huge(0)
end type restart_in_t

! --- Parallel info ---

! Combine information required for non-blocking report loop quantities
! into one type for convenience.
type nb_rep_t
    ! Array to store report loop estimators such as projected energy
    ! etc.
    ! This array must not be deallocated, copied or inspected in any
    ! way in between report loop communication.
    real(dp), allocatable :: rep_info(:)
    ! Array whose entries will contain:
    ! 1. The total number of spawned walkers in a given report loop
    !    which is to be used for calculating the spawning rate.
    ! 2. The number of walkers spawned from a given processor
    !    to all other processors except the current one, which
    !    is used for calculating the total number of walkers for a given
    !    report loop.
    integer :: nb_spawn(2)
    ! Array of requests used for non blocking communications.
    ! This array must not be deallocated, copied or inspected in any
    ! way in between report loop communication. request(nprocs)
    integer, allocatable :: request(:)
end type nb_rep_t

type load_bal_in_t
    ! Number of slots walker lists are initially subdivided into for proc_map
    ! Default = 20. This reverts to 1 when run in serial.
    ! Input option: load_balancing_slots
    integer :: nslots = 20
    ! Population which must be reached before load balancing is attempted.
    ! Default = 1000.
    ! Input option: load_balancing_pop
    integer(int_64) :: pop = 1000
    ! Percentage load imbalance we aim to achieve when performing load balancing.
    ! i.e. min_pop = (1-percent_imbal)*av_pop, max_pop = (1+percent_imbal)*av_pop.
    ! Default = 0.05
    ! Input option: percent_imbal
    real(p) :: percent = 0.05
    ! Maximum number of load balancing attempts.
    ! Default = 2.
    ! Input option: max_load_attempts
    integer :: max_attempts = 2
    ! Write load balancing information every time load balancing is attempted.
    ! Input option: write_load_info
    logical :: write_info = .false.
end type load_bal_in_t

type load_bal_state_t
    ! Tag to check which stage if load balancing is required. This is reset to false
    ! once redistribution of determinants has taken place to ensure load balancing
    ! occurs once during a report loop.
    logical :: needed = .false.
    ! Current number of load balancing attempts.
    integer :: nattempts = 0
    type(proc_map_t) :: proc_map
end type load_bal_state_t

type parallel_t
    ! Type containing information on current state of load imbalance of
    ! particles, including proc_map.
    type(load_bal_state_t) :: load
    ! Type containing arrays necessary for report communication when using
    ! non-blocking communications.
    type(nb_rep_t) :: report_comm
end type parallel_t

!--- Reference determinant ---

type reference_t
    ! Bit string of reference determinant.
    integer(i0), allocatable :: f0(:)
    ! List of occupied orbitals in reference determinant.
    integer, allocatable :: occ_list0(:)
    ! Bit string of reference determinant used to generate the Hilbert space.
    ! This is usually identical to f0, but not necessarily (e.g. if we're doing
    ! a spin-flip calculation).
    integer(i0), allocatable :: hs_f0(:)
    ! hs_f0:hs_occ_list0 as f0:occ_list0.
    integer, allocatable :: hs_occ_list0(:)
    ! CCMC/CIQMC: max number of excitations from the reference to include in
    ! the calculation.
    ! DMQMC: permit density matrix elements to be non-zero only if the two
    ! determinants differ by at most ex_level excitations.
    ! Set to the number of electrons in the system to use the full space.
    integer :: ex_level = -1
    ! Energy of reference determinant.
    real(p) :: H00
    ! Value of <D0|O|D0>, where O is the operator we are sampling.
    ! (Applicable/set only if Hellmann--Feynman sampling is in operation.)
    real(p) :: O00
end type reference_t

! --- semi-stochastic ---

! Types of space.
enum, bind(c)
    ! This option uses an empty deterministic space, and so turns
    ! semi-stochastic off.
    enumerator :: empty_determ_space = 0
    ! This space is generated by choosing the determinants with the highest
    ! populations in the main walker list, at the point the deterministic space
    ! is generated.
    enumerator :: high_pop_determ_space
    ! This option generates the deterministic space by reading the states in
    ! from an HDF5 file (named SEMI.STOCH.*.H5).
    enumerator :: read_determ_space
    ! This option tells the semi-stochastic initialisation routine to use a
    ! deterministic space which has already been created, and which should
    ! be stored in the dets array within the semi_stoch_t instance on input.
    enumerator :: reuse_determ_space
end enum


! Array to hold the indices of deterministic states in the dets array, accessed
! by calculating a hash value. This type is used by the semi_stoch_t type and
! is intended only to be used by this object.
type determ_hash_t
    ! For an example of how to use this type to see if a determinant is in the
    ! deterministic space or not, see the routine check_if_determ.

    !rSeed used in the MurmurHash function to calculate hash values.
    integer :: seed
    ! The size of the hash table (ignoring collisions).
    integer :: nhash
    ! The indicies of the determinants in the semi_stoch_t%dets array.
    ! Note that element nhash+1 should be set equal to determ%tot_size+1.
    ! This helps with avoiding out-of-bounds errors when using this object.
    integer, allocatable :: ind(:) ! (semi_stoch_t%tot_size)
    ! hash_ptr(i) stores the index of the first index in the array ind which
    ! corresponds to a determinant with hash value i.
    ! This is similar to what is done in the CSR sparse matrix type (see
    ! csr.f90).
    integer, allocatable :: hash_ptr(:) ! (nhash+1)
end type determ_hash_t

type semi_stoch_t
    ! True if a semi-stochastic calculation is being performed with this object.
    logical :: doing_semi_stoch = .false.
    ! If true, then the deterministic 'spawning' will be performed in a routine
    ! with an extra MPI call. This routine handles all of the annihilation of
    ! deterministic spawnings with a simple summing of vectors.
    ! If false, then the deterministic spawnings are added to the spawned list
    ! and treated with the standard annihilation routine.
    logical :: separate_annihilation = .false.
    ! Integer to specify which type of deterministic space is being used.
    ! See the various determ_space parameters defined above.
    integer :: space_type = empty_determ_space
    ! The total number of deterministic states on all processes.
    integer :: tot_size = 0
    ! sizes(i) holds the number of deterministic states belonging to process i.
    integer, allocatable :: sizes(:) ! (0:nproc-1)
    ! The Hamiltonian in the deterministic space, stored in a sparse CSR form.
    ! An Hamiltonian element, H_{ij}, is stored in hamil if and only if both
    ! i and j are in the deterministic space.
    type(csrp_t) :: hamil
    ! This array is used to store the values of amplitudes of deterministic
    ! states throughout a QMC calculation.
    real(p), allocatable :: vector(:) ! sizes(iproc)
    ! If separate_annihilation is true, then an extra MPI call is used to join
    ! together the the deterministic vector arrays from each process. This
    ! array is used to hold the results, which will be the list of all
    ! deterministic amplitudes.
    ! If separate_annihilation is not true then this array will remain
    ! deallocated.
    real(p), allocatable :: full_vector(:) ! tot_size
    ! If separate_annihilation is true then this array will hold the indices
    ! of the deterministic states in the main list. This prevents having to
    ! search the whole of the main list for the deterministic states.
    ! If separate_annihilation is not true then this array will remain
    ! deallocated.
    integer, allocatable :: indices(:) ! sizes(iproc)
    ! dets stores the deterministic states across all processes.
    ! All states on process 0 are stored first, then process 1, etc...
    integer(i0), allocatable :: dets(:,:) ! (string_len, tot_size)
    ! A hash table which allows the index of a determinant in dets to be found.
    ! This is done by calculating the hash value of the given determinant.
    type(determ_hash_t) :: hash_table
    ! Deterministic flags of states in the main list. If determ_flags(i) is
    ! equal to 0 then the corresponding state in position i of the main list is
    ! a deterministic state, else it is not.
    integer, allocatable :: flags(:)
    ! Type for holding information about semi-stochastic MPI timings.
    ! This is only used if separate_annihilation is .true.. If it is false
    ! then semi-stochastic communication is performed with the main spawning
    ! communicaton.
    type(parallel_timing_t) :: mpi_time
end type semi_stoch_t

! --- Importance sampling ---

! For the Heisenberg model, several different trial functions can be used in the
! energy estimator. Only a single determinant can be used for the Hubbard model.
enum, bind(c)
    enumerator :: single_basis
    enumerator :: neel_singlet
end enum

! For the Heisenberg model, a guiding function may be used,
! |psi_G> = \sum_{i} a_i |psi_i>, so that the new Hamiltonian matrix elements are
! H_ij^new = (a_i*H_ij)/a_j. This is just importance sampling. These functions
! represent the different types of functions which may be used.
enum, bind(c)
    enumerator :: no_guiding
    ! Note that when we use the Neel singlet state as a guiding function, it must also
    ! be used as the trial function in calculating the projected energy.
    enumerator :: neel_singlet_guiding
end enum

! --- Estimators ---

type particle_t
    ! Current number of walkers stored in the main list (processor dependent).
    ! This is updated during annihilation and merging of the spawned walkers into
    ! the main list.
    integer :: nstates 
    ! Total number of particles on all walkers/determinants (processor dependent)
    ! Updated during death and annihilation and merging.
    ! The first element is the number of normal (Hamiltonian) particles.
    ! Subsequent elements are the number of Hellmann--Feynamnn particles.
    real(p), allocatable :: nparticles(:) ! (sampling_size)
    ! Total number of particles across *all* processors, i.e. \sum_{proc} nparticles_{proc}
    real(p), allocatable :: tot_nparticles(:) ! (sampling_size)
    ! Total number of particles on all determinants for each processor
    real(p), allocatable :: nparticles_proc(:,:) ! (sampling_size,nprocs)
    ! Walker information: main list.
    ! sampling_size is one for each quantity sampled (i.e. 1 for standard
    ! FCIQMC/initiator-FCIQMC, 2 for FCIQMC+Hellmann--Feynman sampling).
    integer :: nspaces
    ! number of additional elements stored for each determinant in dat for
    ! (e.g.) importance sampling.
    integer :: info_size
    ! a) determinants
    integer(i0), allocatable :: states(:,:) ! (string_len, walker_length)
    ! [todo] - nicer referencing for elements of dat and ! pops
    ! b) walker population
    ! NOTE:
    !   When using the real_amplitudes option, pops stores encoded
    !   representations of the true walker populations. To convert
    !   pops(:,i) to the actual population on determinant i, one must
    !   take real(pops(:,i),p)/real_factor. Thus, the resolution
    !   in the true walker populations is 1/real_factor. This is how
    !   non-integers populations are implemented. When not using the real_amplitudes
    !   option, real_factor will be equal to 1, allowing only integer
    !   populations. In general, when one sees that a integer is of kind int_p, it
    !   should be understood that it stores a population in its encoded form.
    integer(int_p), allocatable :: pops(:,:) ! (nspaces,walker_length)
    ! c) Walker information.  This contains:
    ! * Diagonal matrix elements, K_ii.  Storing them avoids recalculation.
    !   K_ii = < D_i | H | D_i > - E_0, where E_0 = <D_0 | H | D_0> and |D_0> is the
    !   reference determinant.  Always the first element.
    ! * Diagonal matrix elements for Hellmann--Feynmann sampling in 2:sampling_size
    !   elements.
    ! * Further data in sampling_size+1:sampling_size:info_size.  For example, when
    !   calculating the projected energy with various trial wavefunctions, it is
    !   useful to store quantites which are expensive to calculate and which are
    !   instead of recalculating them. For the Neel singlet state, the first component
    !   gives the total number of spins up on the first sublattice. The second
    !   component gives the number of 0-1 bonds where the 1 is on the first
    !   sublattice.
    real(p), allocatable :: dat(:,:) ! (sampling_size+info_size,walker_length)
end type particle_t

type spawned_particle_t
    ! List of spawned particles.
    type(spawn_t) :: spawn
    ! List of spawned particles received from other processors.  Used for non-blocking
    ! communications (normal synchronous communications can just use spawn alone).
    type(spawn_t) :: spawn_recv
    ! Rate of spawning.  This is a running total over MC cycles on each processor
    ! until it is summed over processors and averaged over cycles in
    ! update_energy_estimators.
    real(p) :: rspawn
end type spawned_particle_t

type estimators_t
    ! Population of walkers on reference determinant/trace of density matrix.
    real(p) :: D0_population
    ! projected energy
    ! This stores during an FCIQMC report loop
    !   \sum_{i/=0} <D_0|H|D_i> N_i
    ! where D_0 is the reference determinants and N_i is the walker population on
    ! determinant D_i.
    ! The projected energy is given as
    !   <D_0|H|D_0> + \sum_{i/=0} <D_0|H|D_i> N_i/N_0
    ! and so proj_energy must be 'normalised' and averaged over the report loops
    ! accordingly.
    real(p) :: proj_energy
    ! Total number of occupied states across all processors.
    integer :: tot_nstates
    ! The total number of successful spawning events, across all processors.
    integer :: tot_nspawn_events

    ! Hellmann--Feynman sampling (several terms must be accumulated and averaged separately):
    ! Signed population of Hellmann--Feynman particles
    !     \sum_i sign(N_i) \tilde{N}_i,
    ! where
    !     N_i is the Hamiltonian population on |D_i>,
    !     \tilde{N}_i is the Hellmann--Feynman population on |D_i>.
    real(p) :: hf_signed_pop
    ! Population on the reference in the Hellmann--Feynman space, \tilde{N}_0.
    real(p) :: D0_hf_population
    ! \sum_i <D_0|O|D_i> N_i.
    real(p) :: proj_hf_O_hpsip
    ! \sum_i <D_0|H|D_i> \tilde{N}_i.
    real(p) :: proj_hf_H_hfpsip
end type estimators_t

type qmc_state_t
    ! When performing dmqmc calculations, dmqmc_factor = 2.0. This factor is
    ! required because in DMQMC calculations, instead of spawning from one end with
    ! the full probability, we spawn from two different ends with half probability each.
    ! Hence, tau is set to tau/2 in DMQMC calculations, so that an extra factor is not
    ! required in every spawning routine. In the death step however, we use
    ! walker_energies(1,idet), which has a factor of 1/2 included for convenience
    ! already, for conveniece elsewhere. Hence we have to multiply by an extra factor
    ! of 2 to account for the extra 1/2 in tau. dmqmc_factor is set to 1.0 when not
    ! performing a DMQMC calculation, and so can be ignored in these cases.
    real(p) :: dmqmc_factor = 1.0_p
    ! timestep
    real(p) :: tau
    ! number of Monte Carlo cycles done in previous runs (ie from restarts, etc).
    integer :: mc_cycles_done = 0
    ! Energy offset (shift) applied to the Hamiltonian.
    real(p), allocatable :: shift(:) ! (psip_list%nspaces)
    ! The shift is updated at the end of each report loop when vary_shift is true.
    ! When the replica_tricks option is used, the elements
    ! of the shift array refer to the shifts in the corresponding replica systems.
    ! When replica_tricks is not being used, only the first element is used.
    logical, allocatable :: vary_shift(:) ! (psip_list%nspaces)
    ! Convenience handles.
    type(particle_t) :: psip_list
    type(spawned_particle_t) :: spawn_store
    type(reference_t) :: ref
    ! WARNING: par_info is the 'reference/master' (ie correct) version
    ! of parallel_t, in particular of proc_map_t.  However, copies of it
    ! are kept in spawn_t objects, and it is these copies which are used
    ! to determine the processor location of a particle.  It is the programmer's
    ! responsibility to ensure these are kept up to date...
    type(parallel_t) :: par_info
    type(estimators_t) :: estimators
end type qmc_state_t

! Copies of various settings that are required during annihilation.  This avoids having to pass through lots of different
! structs/flags for various settings.  Set in init_qmc.
type annihilation_flags_t
    ! Calculate replicas (ie evolve two wavefunctions/density matrices at once)?
    ! Currently only implemented for DMQMC.
    logical :: replica_tricks = .false.
    ! Propagate a trial density matrix to a specific temeperature.
    logical :: propagate_to_beta = .false.
    ! Trial function used (FCIQMC & Heisenberg model only).
    integer :: trial_function
end type annihilation_flags_t

! --- GLOBAL STATE (TEMPORARY) ---

! Global handle for purity work.
type(qmc_in_t) :: qmc_in_global
type(fciqmc_in_t) :: fciqmc_in_global
type(semi_stoch_in_t) :: semi_stoch_in_global
type(ccmc_in_t) :: ccmc_in_global
type(restart_in_t) :: restart_in_global
type(load_bal_in_t) :: load_bal_in_global

type(annihilation_flags_t) :: annihilation_flags_global

type(particle_t), target :: walker_global

! When using the Neel singlet trial wavefunction, it is convenient
! to store all possible amplitudes in the wavefunction, since
! there are relativley few of them and they are expensive to calculate
real(dp), allocatable :: neel_singlet_amp(:) ! (nsites/2) + 1

! Real amplitudes can be any multiple of 2**(-real_bit_shift). They are
! encoded as integers by multiplying them by 2**(real_bit_shift).
! [todo] - compile-time parameter
integer :: real_bit_shift
! real_factor = 2**(real_bit_shift)
! [todo] - compile-time parameter
integer(int_p) :: real_factor
! [todo] - procedures for encoding and decoding the populations.

end module qmc_data
