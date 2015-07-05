module dmqmc_data

use const, only: p, dp, i0
use hash_table, only: hash_table_t
use spawn_data, only: spawn_t

implicit none

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

! The following are the possible options for the initial density matrices for
! IP-DMQMC.
enum, bind(c)
    ! "Hartree-Fock" density matrix, i.e. \rho = \sum e^{-\beta H_ii} |D_i><D_i|.
    enumerator :: hartree_fock_dm
    ! Free-electron density matrix, i.e. \rho = \sum_i e^{-\beta \sum_j \varepsilon_j \hat{n}_j} |D_i><D_i|.
    enumerator :: free_electron_dm
end enum

! This variable holds the total number of operators which are implemented
! for DMQMC.
integer, parameter :: num_dmqmc_operators = terminator - 1

! This type contains information for a given subsystem. It takes translational
! symmetry into account by storing information for all subsystems which are
! equivalent by translational symmetry.
type subsys_t
    ! The total number of sites in subsystem A.
    integer :: A_nsites
    ! Equivalent to string_len in basis_t, string_len is the length of the byte
    ! array necessary to contain a bit for each subsystem-A basis function.
    integer :: string_len
    ! The sites in subsystem A.
    integer, allocatable :: subsystem_A(:)
    ! B_masks(:,i) has bits set at all bit positions corresponding to sites in
    ! version i of subsystem B, where the different 'versions' correspond to
    ! subsystems which are equivalent by symmetry.
    integer(i0), allocatable :: B_masks(:,:)
    ! bit_pos(i,j,1) contains the position of the bit corresponding to site i in
    ! 'version' j of subsystem A.
    ! bit_pos(i,j,2) contains the element of the bit corresponding to site i in
    ! 'version' j of subsystem A.
    ! Note that site i in a given version is the site that corresponds to site i
    ! in all other versions of subsystem A (and so bit_pos(i,:,1) and
    ! bit_pos(i,:,2) will not be sorted). This is very important for
    ! applications to reduced density matrices, so that equivalent psips will
    ! contribute to the same RDM element.
    integer, allocatable :: bit_pos(:,:,:)
end type subsys_t

type dmqmc_rdm_in_t
    ! The total number of rdms beings calculated (currently only applicable to
    ! instantaneous RDM calculations, not to ground-state RDM calculations,
    ! which only ever calculate one RDM).
    integer :: nrdms = 0

    ! The length of the spawning array for RDMs. Each RDM calculated has the
    ! same length array. Note, this is only used for instantaneous RDMs.
    ! Ground-state RDM calculations allocate an array exactly the size of the
    ! full RDM.
    integer :: spawned_length

    ! If true then the reduced density matricies will be calulated for the 'A'
    ! subsystems specified by the user.
    logical :: doing_rdm = .false.

    ! If true then each subsystem A RDM specified by the user will be accumulated
    ! from the iteration start_averaging until the end of the beat loop, allowing
    ! ground-state estimates of the RDMs to be calculated.
    logical :: calc_ground_rdm = .false.

    ! If true then the reduced density matricies will be calculated for each
    ! subsystem specified by the user at the end of each report loop. These RDMs
    ! can be used to calculate instantaeous estimates at the given beta value.
    ! They are thrown away after these calculation has been performed on them.
    logical :: calc_inst_rdm = .false.

    ! If true then calculate the concurrence for reduced density matrix of two sites.
    logical :: doing_concurrence = .false.

    ! If true then calculate the von Neumann entanglement entropy for specified subsystem.
    logical :: doing_vn_entropy = .false.

    ! If true then the reduced density matrix is output to a file, 'reduced_dm'
    ! each beta loop.
    logical :: output_rdm = .false.

end type dmqmc_rdm_in_t

type dmqmc_in_t
    ! The number of times the program will loop over each value of beta in the main loop.
    integer :: beta_loops = 100

    ! Calculate replicas (ie evolve two wavefunctions/density matrices at once)?
    ! Currently only implemented for DMQMC.
    logical :: replica_tricks = .false.

    ! When performing a ground-state RDM calculation, on what iteration do we
    ! start accumulating the ground-state RDM?
    integer :: start_av_rdm = 0

    ! If this logical is true then the program runs the DMQMC algorithm with
    ! importance sampling.
    ! sampling_probs stores the factors by which the probabilities of
    ! spawning to a larger excitation are reduced by. So, when spawning from
    ! a diagonal element to a element with one excitation, the probability
    ! of spawning is reduced by a factor sampling_probs(1).
    ! accumulated_probs(i) stores the multiplication of all the elements
    ! of sampling_probs up to the ith element. This quantity is often
    ! needed, so it is stored.
    logical :: weighted_sampling = .false.
    ! If vary_weights is true, then instead of using the final sampling
    ! weights for all the iterations, the weights will be gradually increased
    ! until finish_varying_weights, at which point they will be held constant.
    ! altering_factors stores the factors by which each weight is multiplied
    ! at each step.
    logical :: vary_weights = .false.
    ! If this logical is true then the program will calculate the ratios
    ! of the numbers of the psips on neighbouring excitation levels. These
    ! are output so that they can be used when doing importance sampling
    ! for DMQMC, so that each level will have roughly equal numbers of psips.
    ! The resulting new weights are used in the next beta loop.
    logical :: find_weights = .false.
    ! When running a simulation to find some appropriate importance sampling
    ! weights, on which iteration should we start averaging the excitation
    ! distributions (which are used to find these weights)?
    integer :: find_weights_start = 0
    ! If true then the fraction of psips at each
    ! excitation level will be output at each report loop. These fractions
    ! will be stored in the array below.
    ! The number of excitations for a given system is defined by
    ! sys_t%max_number_excitations; see comments in sys_t for more details.
    logical :: calc_excit_dist = .false.
    ! If true then the simulation will start with walkers uniformly distributed
    ! along the diagonal of the entire density matrix, including all symmetry
    ! sectors.
    logical :: all_sym_sectors = .false.
    ! If true then the simulation will start with walkers distributed in all
    ! spin symmetry sectors of the Hamiltonian i.e. all 2S+1 blocks between
    ! -ms...ms.
    logical :: all_spin_sectors = .false.
    ! If half_density_matrix is true then half the density matrix will be
    ! calculated by reflecting spawning onto the lower triangle into the
    ! upper triangle. This is allowed because the density matrix is
    ! symmetric.
    logical :: half_density_matrix = .false.

    ! correlation_sites stores the site positions specified by the users
    ! initially (as orbital labels), to be used in the calculation of
    ! spin correlation functions.
    integer, allocatable :: correlation_sites(:)
    ! correlation_mask is a bit string with a 1 at positions i and j which
    ! are considered when finding the spin correlation function, C(r_{i,j}).
    ! All other bits are set to 0. i and j are chosen by the user initially.
    ! This is not actually an input option but is calculated from
    ! correlation_sites, but we store it here to be pragmatic.
    integer(i0), allocatable :: correlation_mask(:) ! (string_len)

    ! When using the old weighted importance sampling, sampling_probs
    ! stores the factors by which probabilities are to be reduced when spawning
    ! away from the diagonal.
    real(p), allocatable :: sampling_probs(:) ! (max_number_excitations)
    ! When using the old weighted importance sampling, how many iterations are
    ! the weights varied for?
    integer :: finish_varying_weights = 0

    ! Propagate a trial density matrix to a specific temeperature.
    logical :: propagate_to_beta = .false.
    ! Initial density matrix to use in IP-DMQMC see enum at beginning of module
    ! for description of available values.
    integer :: initial_matrix
    ! Use the grand canonical partition function to inititally distribute the psips.
    logical :: grand_canonical_initialisation = .false.
    ! Interpret input init_beta as the inverse reduced temperature, i.e., Beta = 1\Theta = T_F/T.
    logical :: fermi_temperature = .false.
    ! Value of beta which we propagate the density matrix to.
    real(p) :: init_beta = 1.0
    ! Number of metropolis attempts (per psip) we use when generating
    ! the trial density matrix.
    integer :: metropolis_attempts = 0

    ! Input options relating to RDMs in DMQMC.
    type(dmqmc_rdm_in_t) :: rdm

end type dmqmc_in_t

! Spawned lists for rdms.
type rdm_spawn_t
    type(spawn_t) :: spawn
    ! Spawn with the help of a hash table to avoid a sort (which is extremely
    ! expensive when a large number of keys are repeated--seem to hit worst case
    ! performance in quicksort).
    type(hash_table_t) :: ht
end type rdm_spawn_t

!--- Type for all instantaneous RDMs ---
type dmqmc_inst_rdms_t
    ! The total number of rdms beings calculated.
    integer :: nrdms = 0

    ! traces(i,j) holds the trace of replica i of the rdm with label j.
    real(p), allocatable :: traces(:,:) ! (particle_t%nspaces, nrdms)

    type(rdm_spawn_t), allocatable :: spawn(:) ! nrdms

    ! When using the replica_tricks option, if the rdm in the first
    ! simulation if denoted \rho^1 and the ancillary rdm is denoted
    ! \rho^2 then renyi_2 holds:
    ! x = \sum_{ij} \rho^1_{ij} * \rho^2_{ij}.
    ! The indices of renyi_2 hold this value for the various rdms being
    ! calculated. After post-processing averaging, this quantity should
    ! be normalised by the product of the corresponding RDM traces -
    ! call it y. Then the renyi-2 entropy is then given by -log_2(x/y).
    real(p), allocatable :: renyi_2(:) ! (nrdms)
end type dmqmc_inst_rdms_t

!--- Type for a ground state RDM ---
type dmqmc_ground_rdm_t
    ! This stores the reduces matrix, which is slowly accumulated over time
    ! (on each processor).
    real(p), allocatable :: rdm(:,:)
    ! The trace of the ground-state RDM.
    real(p) :: trace
end type dmqmc_ground_rdm_t

type dmqmc_estimates_t
    ! numerators stores the numerators for the estimators in DMQMC. These
    ! are, for a general operator O which we wish to find the thermal average of:
    ! \sum_{i,j} \rho_{ij} * O_{ji}
    ! This variabe will store this value from the first iteration of each
    ! report loop. At the end of a report loop, the values from each
    ! processor are combined and stored in numerators on the parent
    ! processor. This is then output, and the values of numerators
    ! are reset on each processor to start the next report loop.
    real(p) :: numerators(num_dmqmc_operators)

    ! In DMQMC the trace of the density matrix is an important quantity
    ! used in calculating all thermal estimators. This quantity stores
    ! the this value, Tr(\rho), where rho is the density matrix which
    ! the DMQMC algorithm calculates stochastically.
    real(p), allocatable :: trace(:) ! (particle_t%nspaces)

    ! This array is used to hold the number of particles on each excitation
    ! level of the density matrix.
    real(p), allocatable :: excit_dist(:) ! (0:max_number_excitations)

    ! RDM data.

    ! This stores information about the various subsystems corresponding to the
    ! RDMs being studied. Each element of this array corresponds to one of
    ! these RDMs. This is not really an estimate, but we put it here to be
    ! pragmatic.
    type(subsys_t), allocatable :: subsys_info(:) ! (nrdms)

    ! Info about ground-state RDM estimates.
    type(dmqmc_ground_rdm_t) :: ground_rdm

    ! Info about instantaneous temperature-dependent) RDM estimates.
    type(dmqmc_inst_rdms_t) :: inst_rdm
end type dmqmc_estimates_t

!--- Type for weighted sampling parameters ---
type dmqmc_weighted_sampling_t
    ! This holds the factors by which the populations on each excitation level
    ! (from 0 to max_number_excitations) are reduced, relative to DMQMC
    ! without any importance sampling.
    real(p), allocatable :: probs(:) ! (max_number_excitations + 1)
    ! The value of accumulated_probs on the last report cycle.
    real(p), allocatable :: probs_old(:) ! (max_number_excitations + 1)

    ! If varying the weights then this array holds the factors by which the
    ! weights are changed each iteration.
    real(dp), allocatable :: altering_factors(:)
end type dmqmc_weighted_sampling_t

type (dmqmc_in_t) :: dmqmc_in_global
type (dmqmc_estimates_t) :: dmqmc_estimates_global
end module dmqmc_data
