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
    enumerator :: kinetic_ind
    enumerator :: H0_ind
    enumerator :: potential_ind
    enumerator :: HI_ind
    enumerator :: terminator ! unused except in num_dmqmc_operators
   ! NOTE: if you add a new estimator then you must insert it before terminator.
end enum

! The following are the possible options for the initial density matrices for
! IP-DMQMC.
enum, bind(c)
    ! "Hartree-Fock" density matrix, i.e. \rho = \sum e^{-\beta H_ii} |D_i><D_i|.
    enumerator :: hartree_fock_dm = 1
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
    integer :: spawned_length = 0

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
    integer :: initial_matrix = hartree_fock_dm
    ! Use the grand canonical partition function to inititally distribute the psips.
    logical :: grand_canonical_initialisation = .false.
    ! Interpret input target_beta as the inverse reduced temperature, i.e., Beta = 1\Theta = T_F/T.
    logical :: fermi_temperature = .false.
    ! Value of beta which we propagate the density matrix to.
    real(p) :: target_beta = 1.0
    ! Number of metropolis attempts (per psip) we use when generating
    ! the trial density matrix.
    integer :: metropolis_attempts = 0
    ! Do a symmetric version of DMQMC, default true and only changeable for the ip-dmqmc algorithm.
    ! This considerably changes the IP-DMQMC algorithm.
    logical :: symmetric = .true.
    ! Chemical potential used to initialise density matrix.
    real(p) :: chem_pot = 0.0_p

    ! Input options relating to RDMs in DMQMC.
    type(dmqmc_rdm_in_t) :: rdm

    ! Excitation level at which to set a psip to be an initiator.
    ! Can be either set to be negative meaning no initiator space is imposed,
    ! to zero meaning that diagonal elements are automatically initiators or to
    ! two meaning that those density matrix elements at excitation level 2 are also
    ! set to be initiators.
    ! Default: Use normal initiator approximation, i.e., allow initiator space to develop
    ! by itself.
    integer :: initiator_level = -1

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

    ! correlation_mask is a bit string with a 1 at positions i and j which
    ! are considered when finding the spin correlation function, C(r_{i,j}).
    ! All other bits are set to 0. i and j are chosen by the user initially.
    ! Not an estimate, but needed to calculate one.
    integer(i0), allocatable :: correlation_mask(:) ! (string_len)

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
    real(p), allocatable :: probs(:) ! (max_number_excitations + 1 (or 2 if using symmetric ip-dmqmc))
    ! The value of accumulated_probs on the last report cycle.
    real(p), allocatable :: probs_old(:) ! (max_number_excitations + 1 (or 2 if using symmetric ip-dmqmc))
    ! When using the old weighted importance sampling, sampling_probs
    ! stores the factors by which probabilities are to be reduced when spawning
    ! away from the diagonal.
    real(p), allocatable :: sampling_probs(:) ! (max_number_excitations)

    ! If varying the weights then this array holds the factors by which the
    ! weights are changed each iteration.
    real(dp), allocatable :: altering_factors(:)
end type dmqmc_weighted_sampling_t

contains


    subroutine dmqmc_in_t_json(js, dmqmc, terminal)

        ! Serialise a dmqmc_in_t object in JSON format.

        ! In/Out:
        !   js: json_out_t controlling the output unit and handling JSON internal state.  Unchanged on output.
        ! In:
        !   dmqmc_in: dmqmc_in_t object containing dmqmc input values (including any defaults set).
        !   terminal (optional): if true, this is the last entry in the enclosing JSON object.  Default: false.

        use json_out

        type(json_out_t), intent(inout) :: js
        type(dmqmc_in_t), intent(in) :: dmqmc
        logical, intent(in), optional :: terminal

        call json_object_init(js, 'dmqmc')
        call json_write_key(js, 'beta_loops', dmqmc%beta_loops)
        call json_write_key(js, 'replica_tricks', dmqmc%replica_tricks)
        call json_write_key(js, 'start_av_rdm', dmqmc%start_av_rdm)
        call json_write_key(js, 'weighted_sampling', dmqmc%weighted_sampling)
        call json_write_key(js, 'vary_weights', dmqmc%vary_weights)
        call json_write_key(js, 'find_weights', dmqmc%find_weights)
        call json_write_key(js, 'find_weights_start', dmqmc%find_weights_start)
        call json_write_key(js, 'calc_excit_dist', dmqmc%calc_excit_dist)
        call json_write_key(js, 'all_sym_sectors', dmqmc%all_sym_sectors)
        call json_write_key(js, 'all_spin_sectors', dmqmc%all_spin_sectors)
        call json_write_key(js, 'initiator_level', dmqmc%initiator_level)
        if (allocated(dmqmc%sampling_probs)) then
            call json_write_key(js, 'sampling_probs', dmqmc%sampling_probs)
        else
            call json_write_key(js, 'sampling_probs', '[]')
        end if
        call json_write_key(js, 'finish_varying_weights', dmqmc%finish_varying_weights)
        call json_write_key(js, 'fermi_temperature', dmqmc%fermi_temperature)
        call json_write_key(js, 'target_beta', dmqmc%target_beta, terminal=.true.)
        call json_object_end(js, terminal)

    end subroutine dmqmc_in_t_json

    subroutine ipdmqmc_in_t_json(js, dmqmc, terminal)

        ! Serialise a dmqmc_in_t object in JSON format. IP-DMQMC specific input
        ! options.

        ! In/Out:
        !   js: json_out_t controlling the output unit and handling JSON internal state.  Unchanged on output.
        ! In:
        !   dmqmc_in: dmqmc_in_t object containing dmqmc input values (including any defaults set).
        !   terminal (optional): if true, this is the last entry in the enclosing JSON object.  Default: false.

        use json_out

        type(json_out_t), intent(inout) :: js
        type(dmqmc_in_t), intent(in) :: dmqmc
        logical, intent(in), optional :: terminal

        call json_object_init(js, 'ipdmqmc')
        call json_write_key(js, 'propagate_to_beta', dmqmc%propagate_to_beta)
        select case(dmqmc%initial_matrix)
        case(hartree_fock_dm)
            call json_write_key(js, 'initial_matrix', 'hartree_fock')
        case(free_electron_dm)
            call json_write_key(js, 'initial_matrix', 'free_electron')
        case default
            call json_write_key(js, 'initial_matrix', dmqmc%initial_matrix)
        end select
        call json_write_key(js, 'grand_canonical_initialisation', dmqmc%grand_canonical_initialisation)
        call json_write_key(js, 'symmetric', dmqmc%symmetric)
        call json_write_key(js, 'chem_pot', dmqmc%chem_pot)
        call json_write_key(js, 'metropolis_attempts', dmqmc%metropolis_attempts, terminal=.true.)
        call json_object_end(js, terminal)

    end subroutine ipdmqmc_in_t_json

    subroutine rdm_in_t_json(js, rdm, terminal)

        ! Serialise a dmqmc_rdm_in_t object in JSON format.
        ! options.

        ! In/Out:
        !   js: json_out_t controlling the output unit and handling JSON internal state.  Unchanged on output.
        ! In:
        !   rdm_in: dmqmc_rdm_in_t object containing dmqmc input values (including any defaults set).
        !   terminal (optional): if true, this is the last entry in the enclosing JSON object.  Default: false.

        use json_out

        type(json_out_t), intent(inout) :: js
        type(dmqmc_rdm_in_t), intent(in) :: rdm
        logical, intent(in), optional :: terminal

        call json_object_init(js, 'rdm')
        call json_write_key(js, 'nrdms', rdm%nrdms)
        call json_write_key(js, 'spawned_length', rdm%spawned_length)
        call json_write_key(js, 'doing_rdm', rdm%doing_rdm)
        call json_write_key(js, 'calc_ground_rdm', rdm%calc_ground_rdm)
        call json_write_key(js, 'calc_inst_rdm', rdm%calc_inst_rdm)
        call json_write_key(js, 'doing_concurrence', rdm%doing_concurrence)
        call json_write_key(js, 'doing_vn_entropy', rdm%doing_vn_entropy)
        call json_write_key(js, 'output_rdm', rdm%output_rdm, terminal=.true.)
        call json_object_end(js, terminal)

    end subroutine rdm_in_t_json

    subroutine operators_in_t_json(js, terminal)

        ! serialise operators in json format.
        ! options.

        ! in/out:
        !   js: json_out_t controlling the output unit and handling json internal state.  unchanged on output.
        ! in:
        !   terminal (optional): if true, this is the last entry in the enclosing json object.  default: false.

        use json_out
        use calc, only: doing_dmqmc_calc, dmqmc_energy, dmqmc_energy_squared, dmqmc_correlation, &
                        dmqmc_staggered_magnetisation, dmqmc_rdm_r2, dmqmc_full_r2, dmqmc_kinetic_energy, &
                        dmqmc_potential_energy, dmqmc_H0_energy, dmqmc_HI_energy

        type(json_out_t), intent(inout) :: js
        logical, intent(in), optional :: terminal

        call json_object_init(js, 'operators')
        call json_write_key(js, 'energy', doing_dmqmc_calc(dmqmc_energy))
        call json_write_key(js, 'energy_squared', doing_dmqmc_calc(dmqmc_energy_squared))
        call json_write_key(js, 'kinetic_energy', doing_dmqmc_calc(dmqmc_kinetic_energy))
        call json_write_key(js, 'potential_energy', doing_dmqmc_calc(dmqmc_potential_energy))
        call json_write_key(js, 'H0_energy', doing_dmqmc_calc(dmqmc_H0_energy))
        call json_write_key(js, 'HI_energy', doing_dmqmc_calc(dmqmc_HI_energy))
        call json_write_key(js, 'correlation_fn', doing_dmqmc_calc(dmqmc_correlation))
        call json_write_key(js, 'staggered_mad_ind', doing_dmqmc_calc(dmqmc_staggered_magnetisation))
        call json_write_key(js, 'rdm_r2', doing_dmqmc_calc(dmqmc_rdm_r2))
        call json_write_key(js, 'full_r2', doing_dmqmc_calc(dmqmc_full_r2), terminal=.true.)
        call json_object_end(js, terminal)

    end subroutine operators_in_t_json

    subroutine subsys_t_json(js, subsys, terminal)

        ! Serialise an array of subsys_t objects in JSON format.

        ! In/Out:
        !   js: json_out_t controlling the output unit and handling JSON internal state.  Unchanged on output.
        ! In:
        !   subsys: array of subsys_t objects.  The subsystem_A array is printed out for each subsys_t object,
        !       with each object labelled by its index in the "subsys" JSON object.
        !   terminal (optional): if true, this is the last entry in the enclosing JSON object.  Default: false.

        use json_out
        use utils, only: int_fmt

        type(json_out_t), intent(inout) :: js
        type(subsys_t), intent(in) :: subsys(:)
        logical, intent(in), optional :: terminal
        integer :: i
        character(10) :: ic

        ! Only need to output subsystem_A to reproduce the subsystem information.
        call json_object_init(js, 'subsys')
        do i = 1, size(subsys)
            write (ic,"("//int_fmt(i, 0)//")") i
            call json_write_key(js, trim(ic), subsys(i)%subsystem_A, i==size(subsys))
        end do
        call json_object_end(js, terminal)

    end subroutine subsys_t_json

    subroutine write_dmqmc_report_header(ntypes, dmqmc_in, max_excit)

        ! Write header for DMQMC specific information.

        ! In:
        !    ntypes: number of particle types being sampled.
        !    dmqmc_in: input options for dmqmc calculations.
        !    max_excit: maximum number of excitations in system.

        use calc, only: doing_calc, hfs_fciqmc_calc, dmqmc_calc, doing_dmqmc_calc
        use calc, only: dmqmc_energy, dmqmc_energy_squared, dmqmc_staggered_magnetisation
        use calc, only: dmqmc_correlation, dmqmc_full_r2, dmqmc_rdm_r2, dmqmc_kinetic_energy
        use calc, only: dmqmc_H0_energy, dmqmc_potential_energy, dmqmc_HI_energy
        use utils, only: int_fmt

        integer, intent(in) :: ntypes
        type(dmqmc_in_t), intent(in) :: dmqmc_in
        integer, intent(in) :: max_excit

        integer :: i, j
        character(16) :: excit_header

        write (6,'(1X,"Information printed out every QMC report loop:",/)')
        write (6,'(1X,"Shift: the energy offset calculated at the end of the report loop.")')
        if (doing_dmqmc_calc(dmqmc_full_r2)) then
            write (6, '(1X,a104)') 'Trace: The current total population on the diagonal elements of the &
                                 &first replica of the density matrix.'
            write (6, '(1X,a107)') 'Trace 2: The current total population on the diagonal elements of the &
                                 &second replica of the density matrix.'
        else
            write (6, '(1X,a83)') 'Trace: The current total population on the diagonal elements of the &
                                 &density matrix.'
        end if
        if (doing_dmqmc_calc(dmqmc_full_r2)) then
            write (6, '(1X,a81)') 'Full S2: The numerator of the estimator for the Renyi entropy of the &
                                  &full system.'
        end if
        if (doing_dmqmc_calc(dmqmc_energy)) then
            write (6, '(1X,a92)') '\sum\rho_{ij}H_{ji}: The numerator of the estimator for the expectation &
                                 &value of the energy.'
        end if
        if (doing_dmqmc_calc(dmqmc_energy_squared)) then
            write (6, '(1X,a100)') '\sum\rho_{ij}H2{ji}: The numerator of the estimator for the expectation &
                                 &value of the energy squared.'
        end if
        if (doing_dmqmc_calc(dmqmc_correlation)) then
            write (6, '(1X,a111)') '\sum\rho_{ij}S_{ji}: The numerator of the estimator for the expectation &
                                 &value of the spin correlation function.'
        end if
        if (doing_dmqmc_calc(dmqmc_staggered_magnetisation)) then
            write (6, '(1X,a109)') '\sum\rho_{ij}M2{ji}: The numerator of the estimator for the expectation &
                                 &value of the staggered magnetisation.'
        end if
        if (doing_dmqmc_calc(dmqmc_rdm_r2)) then
            write (6, '(1x,a73)') 'RDM(n) S2: The numerator of the estimator for the Renyi entropy of RDM n.'
        end if
        if (dmqmc_in%rdm%calc_inst_rdm) then
            write (6, '(1x,a83)') 'RDM(n) trace m: The current total population on the diagonal of replica m &
                                  &of RDM n.'
        end if
        if (dmqmc_in%calc_excit_dist) write (6, '(1x,a86)') &
                'Excit. level n: The fraction of particles on excitation level n of the density matrix.'

        write (6,'(1X,"# particles: current total population of Hamiltonian particles.")')
        write (6,'(1X,"# states: number of many-particle states occupied.")')
        write (6,'(1X,"# spawn_events: number of successful spawning events across all processors.")')
        write (6,'(1X,"R_spawn: average rate of spawning across all processors.")')
        write (6,'(1X,"time: average time per Monte Carlo cycle.",/)')
        write (6,'(1X,"Note that all particle populations are averaged over the report loop.",/)')

        ! Header of data table.
        write (6,'(1X,a12,3X,a13,17X,a5)', advance = 'no') '# iterations','Instant shift','Trace'
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
        if (doing_dmqmc_calc(dmqmc_kinetic_energy)) then
            write (6, '(2X,a19)', advance = 'no') '\sum\rho_{ij}T_{ji}'
        end if
        if (doing_dmqmc_calc(dmqmc_H0_energy)) then
            write (6, '(2X,a19)', advance = 'no') '\sum\rho_{ij}H0{ji}'
        end if
        if (doing_dmqmc_calc(dmqmc_HI_energy)) then
            write (6, '(2X,a19)', advance = 'no') '\sum\rho_{ij}HI{ji}'
        end if
        if (doing_dmqmc_calc(dmqmc_potential_energy)) then
            write (6, '(2X,a19)', advance = 'no') '\sum\rho_{ij}U_{ji}'
        end if
        if (doing_dmqmc_calc(dmqmc_rdm_r2)) then
            do i = 1, dmqmc_in%rdm%nrdms
                write (6, '(16X,a3,'//int_fmt(i,0)//',1x,a2)', advance = 'no') 'RDM', i, 'S2'
            end do
        end if
        if (dmqmc_in%rdm%calc_inst_rdm) then
            do i = 1, dmqmc_in%rdm%nrdms
                do j = 1, ntypes
                    write (6, '(7X,a3,'//int_fmt(i,0)//',1x,a5,1x,'//int_fmt(j,0)//')', advance = 'no') &
                            'RDM', i, 'trace', j
                end do
            end do
        end if
        if (dmqmc_in%calc_excit_dist) then
            do i = 0, max_excit
                write (excit_header, '("Excit. level",1X,'//int_fmt(i,0)//')') i
                write (6, '(5X,a16)', advance='no') excit_header
            end do
        end if

        write (6, '(3X,a11,6X)', advance='no') '# particles'
        write (6,'(3X,"# states  # spawn_events  R_spawn    time")')

    end subroutine write_dmqmc_report_header

    subroutine write_dmqmc_report(qmc_in, qs, ireport, ntot_particles, elapsed_time, comment, &
                                  dmqmc_in, dmqmc_estimates)

        ! Write the report line at the end of a report loop.

        ! In:
        !    qmc_in: input options relating to QMC methods.
        !    qs: QMC state (containing shift and various estimators).
        !    ireport: index of the report loop.
        !    ntot_particles: total number of particles in main walker list.
        !    elapsed_time: time taken for the report loop.
        !    comment: if true, then prefix the line with a #.
        !    dmqmc_in: input options relating to DMQMC.
        !    dmqmc_estimates: type containing all DMQMC estimates to be printed.

        use calc, only: doing_calc, dmqmc_calc, doing_dmqmc_calc
        use calc, only: dmqmc_energy, dmqmc_energy_squared, dmqmc_full_r2, dmqmc_rdm_r2
        use calc, only: dmqmc_correlation, dmqmc_staggered_magnetisation, dmqmc_kinetic_energy
        use calc, only: dmqmc_H0_energy, dmqmc_potential_energy, dmqmc_HI_energy
        use qmc_data, only: qmc_in_t, qmc_state_t

        type(qmc_in_t), intent(in) :: qmc_in
        type(qmc_state_t), intent(in) :: qs
        integer, intent(in) :: ireport
        real(dp), intent(in) :: ntot_particles(:)
        real, intent(in) :: elapsed_time
        logical, intent(in) :: comment
        type(dmqmc_in_t), intent(in) :: dmqmc_in
        type(dmqmc_estimates_t), intent(in) :: dmqmc_estimates

        integer :: mc_cycles, i, j, ntypes

        ntypes = size(ntot_particles)

        mc_cycles = ireport*qmc_in%ncycles

        if (comment) then
            write (6,'(1X,"#",1X)', advance='no')
        else
            write (6,'(3X)', advance='no')
        end if

        write (6,'(i10,2X,es17.10,2X,es17.10)',advance = 'no') &
            (qs%mc_cycles_done+mc_cycles-qmc_in%ncycles), qs%shift(1), dmqmc_estimates%trace(1)
        ! The trace on the second replica.
        if (doing_dmqmc_calc(dmqmc_full_r2)) then
            write(6, '(3X,es17.10)',advance = 'no') dmqmc_estimates%trace(2)
        end if

        ! Renyi-2 entropy for the full density matrix.
        if (doing_dmqmc_calc(dmqmc_full_r2)) then
            write (6, '(4X,es17.10)', advance = 'no') dmqmc_estimates%numerators(full_r2_ind)
        end if

        ! Energy.
        if (doing_dmqmc_calc(dmqmc_energy)) then
            write (6, '(4X,es17.10)', advance = 'no') dmqmc_estimates%numerators(energy_ind)
        end if

        ! Energy squared.
        if (doing_dmqmc_calc(dmqmc_energy_squared)) then
            write (6, '(4X,es17.10)', advance = 'no') dmqmc_estimates%numerators(energy_squared_ind)
        end if

        ! Correlation function.
        if (doing_dmqmc_calc(dmqmc_correlation)) then
            write (6, '(4X,es17.10)', advance = 'no') dmqmc_estimates%numerators(correlation_fn_ind)
        end if

        ! Staggered magnetisation.
        if (doing_dmqmc_calc(dmqmc_staggered_magnetisation)) then
            write (6, '(4X,es17.10)', advance = 'no') dmqmc_estimates%numerators(staggered_mag_ind)
        end if

        ! Kinetic energy
        if (doing_dmqmc_calc(dmqmc_kinetic_energy)) then
            write (6, '(4X,es17.10)', advance = 'no') dmqmc_estimates%numerators(kinetic_ind)
        end if

        ! H^0 energy, where H = H^0 + V.
        if (doing_dmqmc_calc(dmqmc_H0_energy)) then
            write (6, '(4X,es17.10)', advance = 'no') dmqmc_estimates%numerators(H0_ind)
        end if

        ! H^I energy, where H^I = exp(-(beta-tau)/2 H^0) H exp(-(beta-tau)/2. H^0).
        if (doing_dmqmc_calc(dmqmc_HI_energy)) then
            write (6, '(4X,es17.10)', advance = 'no') dmqmc_estimates%numerators(HI_ind)
        end if

        ! Potential energy.
        if (doing_dmqmc_calc(dmqmc_potential_energy)) then
            write (6, '(4X,es17.10)', advance = 'no') dmqmc_estimates%numerators(potential_ind)
        end if

        ! Renyi-2 entropy for all RDMs being sampled.
        if (doing_dmqmc_calc(dmqmc_rdm_r2)) then
            do i = 1, dmqmc_in%rdm%nrdms
                write (6, '(6X,es17.10)', advance = 'no') dmqmc_estimates%inst_rdm%renyi_2(i)
            end do
        end if

        ! Traces for instantaneous RDM estimates.
        if (dmqmc_in%rdm%calc_inst_rdm) then
            do i = 1, dmqmc_in%rdm%nrdms
                do j = 1, ntypes
                    write (6, '(2x,es17.10)', advance = 'no') dmqmc_estimates%inst_rdm%traces(j,i)
                end do
            end do
        end if

        ! The distribution of walkers on different excitation levels of the
        ! density matrix.
        if (dmqmc_in%calc_excit_dist) then
            do i = 0, ubound(dmqmc_estimates%excit_dist,1)
                write (6, '(4X,es17.10)', advance = 'no') dmqmc_estimates%excit_dist(i)/ntot_particles(1)
            end do
        end if

        write (6, '(2X,es17.10)', advance='no') ntot_particles(1)
        write (6,'(2X,i10,4X,i12,2X,f7.4,2X,f7.3)') qs%estimators%tot_nstates, qs%estimators%tot_nspawn_events, &
                                             qs%spawn_store%rspawn, elapsed_time/qmc_in%ncycles

    end subroutine write_dmqmc_report

end module dmqmc_data
