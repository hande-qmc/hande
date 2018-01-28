module excit_gens
use const
implicit none

private
public :: alloc_alias_table_data_t, dealloc_excit_gen_data_t, p_single_double_t, excit_gen_power_pitzer_t, excit_gen_heat_bath_t
public :: move_pattempt_data, excit_gen_data_t

!Data for the power_pitzer/heat bath excit gens

! The integer types have been chosen to be int_32 as they never need to index more than 2^31-1 basis functions.
! WARNING: If int_bas gets modified, MPI calls have to be modified as well (in other files).
integer(int_32), parameter :: int_bas = int_32

! [review] - AJWT: What is this type for?
type p_single_double_coll_t
    real(p) :: h_pgen_singles_sum = 0.0_p ! hmatel/pgen sum for singles
    real(p) :: h_pgen_doubles_sum = 0.0_p ! hamtel/pgen sum for doubles
    real(p) :: excit_gen_singles = 0.0_p ! number of valid singles excitations created
    real(p) :: excit_gen_doubles = 0.0_p ! number of valid doubles excitations created
end type p_single_double_coll_t

type p_single_double_t
    ! collects data to update pattempt_single (and therefore pattempt_double) such that
    ! the means in the distribution of hmatel/pgen for singles and doubles are roughly equal.
    type(p_single_double_coll_t) :: total ! the same on all MPI procs, stored in restart file.
![review] - AJWT: tmp is quite a bad name.  accum?  accum_proc?
    type(p_single_double_coll_t) :: tmp ! gets zeroed after each report loop, different on each MPI proc. Gets added onto total.
    ! [todo] - find a way to make the next two variables parameters inside types.
    ! [todo] - Be careful when changing them - if restarting from a legacy file that might be confusing.
    real(p) :: every_attempts = 10000.0_p ! update pattempt_single every every_attempts a single or double ex. happened
    real(p) :: every_min_attempts = 10.0_p ! unless there were not more than every_min_attempts of single or double ex. since
    real(p) :: counter = 1.0_p
    logical :: vary_psingles = .false. ! still update pattempt_singles
    logical :: pattempt_restart_store = .false. ! has pattempt_* information from p_single_double been stored when restarting
    logical :: overflow_loc = .false. ! Does adding a 1 still increase the number (excit_gen_singles and excit_gen_doubles)?
end type p_single_double_t

type alias_table_data_one_ind_t
    ! aliasU(:) stores the alias table of real U for selecting an orbital.
    real(p), allocatable :: aliasU(:) 
    ! aliasK(:) stores the alias table of integer K for selecting an orbital.
    integer(int_bas), allocatable :: aliasK(:) 
    ! weights(:) stores the weights for selecting an orbital.
    real(p), allocatable :: weights(:) 
    ! weights_tot is the sum over j in weights(j). 
    real(p)  :: weights_tot
end type alias_table_data_one_ind_t

type alias_table_data_two_ind_t
    ! e.g. aliasU(:,i) stores the alias table of real U considering an excitation from i 
    ! to another orbital. 
    real(p), allocatable :: aliasU(:,:) 
    ! aliasK(:,i) stores the alias table of integer K considering an excitation from
    ! i to another orbital.
    integer(int_bas), allocatable :: aliasK(:,:) 
    ! weights(:,i) stores the weights considering an excitation from i to another orbital.
    real(p), allocatable :: weights(:,:) 
    ! weights_tot(i) is the sum over j in weights(j,i). 
    real(p), allocatable :: weights_tot(:) 
end type alias_table_data_two_ind_t

type alias_table_data_three_ind_t
    ! e.g. aliasU(:,syma,i) stores the alias table of real U considering an excitation
    ! from i to any orbital of symmetry syma.
    real(p), allocatable :: aliasU(:,:,:) 
    ! aliasK(:,syma,i) stores the alias table of integer K considering an excitation
    ! from i to any orbital of symmetry syma.
    integer(int_bas), allocatable :: aliasK(:,:,:) 
    ! weights(:,syma,i) stores the weights considering an excitation from i to
    ! orbital of symmetry syma.
    real(p), allocatable :: weights(:,:,:) 
    ! weights_tot(syma,i) is the sum over j in weights(j,syma,i).
    real(p), allocatable :: weights_tot(:,:)
end type alias_table_data_three_ind_t

type alias_table_data_four_ind_t
    ! aliasU(:,a,j,i) stores the alias table of real U considering an excitation
    ! to an orbital having selected ija.
    real(p), allocatable :: aliasU(:,:,:,:) 
    ! aliasK(:,syma,i) stores the alias table of integer K considering an excitation
    ! to an orbital having selected ija.
    integer(int_bas), allocatable :: aliasK(:,:,:,:) 
    ! weights(:,syma,i) stores the weights considering an excitation from i to
    ! an orbital having selected ija.
    real(p), allocatable :: weights(:,:,:,:)
    ! weights_tot(a,j,i) is the sum over b in weights(b,a,j,i).
    real(p), allocatable :: weights_tot(:,:,:)
end type alias_table_data_four_ind_t

! Type containing alias tables, etc, needed when using the Power Pitzer ("occ_ref") and the
! Power Pitzer Order N excitation generators.
type excit_gen_power_pitzer_t
    ! Minimum ratio over number of orbitals of individual weight/total weight. Aims to reduce p_gens that are too big.
    real(p) :: power_pitzer_min_weight
    ! Alias table information for the Power Pitzer excitation generator, selecting a from i in a double excitation.
    type(alias_table_data_two_ind_t) :: pp_ia_d
    ! Alias table information for the Power Pitzer excitation generator, selecting b from j in a double excitation.
    type(alias_table_data_three_ind_t) :: pp_jb_d
    ! Alias table information for the Power Pitzer order N excitation generator, selecting i in a double excitation.
    type(alias_table_data_one_ind_t) :: ppn_i_d
    ! Alias table information for the Power Pitzer order N excitation generator, selecting i in a single excitation.
    type(alias_table_data_one_ind_t) :: ppn_i_s
    ! Alias table information for the Power Pitzer order N excitation generator, selecting j (after i) in a double excitation.
    type(alias_table_data_two_ind_t) :: ppn_ij_d
    ! Alias table information for the Power Pitzer order N excitation generator, selecting a from i in a double excitation.
    type(alias_table_data_two_ind_t) :: ppn_ia_d
    ! Alias table information for the Power Pitzer order N excitation generator, selecting a from i in a single excitation.
    type(alias_table_data_two_ind_t) :: ppn_ia_s
    ! Alias table information for the Power Pitzer order N excitation generator, selecting b from j in a double excitation.
    type(alias_table_data_three_ind_t) :: ppn_jb_d
    ! virt_list_alpha is a list of all virtual alpha orbitals as seen from the reference.
    ! Length of array: (sys%nvirt_alpha)
    integer(int_bas), allocatable :: virt_list_alpha(:) 
    ! virt_list_beta is a list of all virtual beta orbitals as seen from the reference.
    ! Length of array: (sys%nvirt_beta)
    integer(int_bas), allocatable :: virt_list_beta(:) 
    ! occ_list(:) is the list of occupied orbitals in the reference.
    ! Length of array: (sys%nel+1) - The +1 is a pad.
    integer(int_bas), allocatable :: occ_list(:)
    ! List of all alpha orbitals.
    ! Length of array: sys%basis%nbasis 
    integer(int_bas), allocatable :: all_list_alpha(:)
    ! List of all beta orbitals.
    ! Length of array: sys%basis%nbasis
    integer(int_bas), allocatable :: all_list_beta(:)
    ! [todo] - put in sys?
    ! Number of alpha spinorbitals
    integer :: n_all_alpha
    ! Number of beta spinorbitals
    integer :: n_all_beta
end type excit_gen_power_pitzer_t

! Type containing alias tables, etc, needed when using the heat bath excitation generators.
type excit_gen_heat_bath_t
    ! Length of array: sys%basis%nbasis
    real(p), allocatable :: i_weights(:) ! This will hold S_p in the Holmes JCTC paper.
    ! Length of array: sys%basis%nbasis, sys%basis%nbasis
    real(p), allocatable :: ij_weights(:,:) ! This stores D_pq.
    ! Alias table information for selecting a from ij in a double excitation/when using the original
    ! implementation.
    type(alias_table_data_three_ind_t) :: hb_ija
    ! Alias table information for selecting b from ija in a double excitation.
    type(alias_table_data_four_ind_t) :: hb_ijab
end type excit_gen_heat_bath_t

!Type containing data for excitation generators
type excit_gen_data_t
    ! Excitation generator to use, duplicated from qmc_in.
    integer :: excit_gen

    ! Probability of attempting single or double excitations.
    real(p) :: pattempt_single, pattempt_double
    ! Probability of selecting ij to be parallel. Used in no_renorm_spin, renorm_spin
    ! excitation generators.
    real(p) :: pattempt_parallel

    ! When creating an arbitrary excitation, k_i,k_j->k_a,k_b, we must conserve
    ! crystal momentum, k_i+k_j-k_a-k_b=0.  Hence once we've chosen k_i, k_j and
    ! k_a, k_b is uniquely defined.  Further, once we've chosen k_i and k_j and if
    ! we require k_b to exist in the basis, then only certain values of k_a are
    ! permitted.  sys%ueg%ternary_conserve(0,k1,k2,k3) gives how many k_a are permitted
    ! for k_i+k_j = (k1,k2,k3) and sys%ueg%ternary_conserve(1:,k1,k2,k3) gives a bit
    ! string with only bytes set corresponding to permitted k_a values.  Note only
    ! basis functions corresponding to *alpha* orbitals are set.
    ! For systems with dimensionality lower than 3, the higher ki values are set to
    ! 0, i.e. dimensions:
    ! (0:bit_string_len,-N:N,0,0) (1D)
    ! (0:bit_string_len,-N:N,-N:N,0) (2D)
    ! (0:bit_string_len,-N:N,-N:N,-N:N) (3D)
    ! NOTE: this contains values of k_i+k_j which cannot be formed by the basis with
    ! the energy cutoff.  Memory can be saved by not using a cubic array for
    ! k_i+k_j...
    integer(i0), allocatable :: ueg_ternary_conserve(:,:,:,:)
    type(excit_gen_power_pitzer_t) :: excit_gen_pp
    type(excit_gen_heat_bath_t) :: excit_gen_hb
    type (p_single_double_t) :: p_single_double
end type excit_gen_data_t

interface alloc_alias_table_data_t
    module procedure :: alloc_alias_table_data_one_ind_t
    module procedure :: alloc_alias_table_data_two_ind_t
    module procedure :: alloc_alias_table_data_three_ind_t
    module procedure :: alloc_alias_table_data_2indarray_three_ind_t
    module procedure :: alloc_alias_table_data_four_ind_t
end interface

interface dealloc_alias_table_data_t
    module procedure :: dealloc_alias_table_data_one_ind_t
    module procedure :: dealloc_alias_table_data_two_ind_t
    module procedure :: dealloc_alias_table_data_three_ind_t
    module procedure :: dealloc_alias_table_data_four_ind_t
end interface

contains
    
    subroutine move_pattempt_data(excit_gen_data_old, excit_gen_data_new)

        ! [todo] - Move more excit_gen_data?
        ! Move allocated memory from one excit_gen_data_t to another.

        ! In/Out:
        !   excit_gen_data_old: assigned excit_gen_data_t;
        ! Out:
        !   excit_gen_data_new: newly assigned excit_gen_data_t;

        type(excit_gen_data_t), intent(inout) :: excit_gen_data_old
        type(excit_gen_data_t), intent(out) :: excit_gen_data_new
        
        excit_gen_data_new%pattempt_single = excit_gen_data_old%pattempt_single
        excit_gen_data_new%pattempt_double = excit_gen_data_old%pattempt_double
        excit_gen_data_new%p_single_double%counter = excit_gen_data_old%p_single_double%counter
        ! [todo] - Below a type gets set equal. Check.
        excit_gen_data_new%p_single_double%total = excit_gen_data_old%p_single_double%total
        
    end subroutine move_pattempt_data

    subroutine dealloc_excit_gen_data_t(excit_gen_data)

        ! Deallocate the excitation generator data.

        ! In/Out:
        !   excit_gen_data: excit_gen_data_t to be deallocated.

        ! It is a bit superfluous to have this subroutine for one deallocation, but it will be more
        ! useful in future.

        type(excit_gen_data_t), intent(inout) :: excit_gen_data

        if (allocated(excit_gen_data%ueg_ternary_conserve)) deallocate(excit_gen_data%ueg_ternary_conserve)

        call dealloc_excit_gen_power_pitzer_t(excit_gen_data%excit_gen_pp)
        call dealloc_excit_gen_heat_bath_t(excit_gen_data%excit_gen_hb)

    end subroutine dealloc_excit_gen_data_t

    subroutine dealloc_excit_gen_power_pitzer_t(excit_gen_pp)
        
        ! Deallocate the power pitzer excitation generator data.

        ! In/Out:
        !   excit_gen_pp: excit_gen_power_pitzer_t to be deallocated.
        use checking, only: check_deallocate

        type(excit_gen_power_pitzer_t), intent(inout) :: excit_gen_pp
        integer :: ierr

        call dealloc_alias_table_data_t(excit_gen_pp%pp_ia_d)
        call dealloc_alias_table_data_t(excit_gen_pp%pp_jb_d)
        call dealloc_alias_table_data_t(excit_gen_pp%ppn_i_d)
        call dealloc_alias_table_data_t(excit_gen_pp%ppn_i_s)
        call dealloc_alias_table_data_t(excit_gen_pp%ppn_ij_d)
        call dealloc_alias_table_data_t(excit_gen_pp%ppn_ia_d)
        call dealloc_alias_table_data_t(excit_gen_pp%ppn_ia_s)
        call dealloc_alias_table_data_t(excit_gen_pp%ppn_jb_d)
        
        if (allocated(excit_gen_pp%virt_list_alpha)) then
            deallocate(excit_gen_pp%virt_list_alpha, stat=ierr)
            call check_deallocate('excit_gen_pp%virt_list_alpha',ierr)
        end if
        if (allocated(excit_gen_pp%virt_list_beta)) then
            deallocate(excit_gen_pp%virt_list_beta, stat=ierr)
            call check_deallocate('excit_gen_pp%virt_list_beta',ierr)
        end if
        if (allocated(excit_gen_pp%occ_list)) then
            deallocate(excit_gen_pp%occ_list, stat=ierr)
            call check_deallocate('excit_gen_pp%occ_list',ierr)
        end if
        if (allocated(excit_gen_pp%all_list_alpha)) then
            deallocate(excit_gen_pp%all_list_alpha, stat=ierr)
            call check_deallocate('excit_gen_pp%all_list_alpha',ierr)
        end if
        if (allocated(excit_gen_pp%all_list_beta)) then
            deallocate(excit_gen_pp%all_list_beta, stat=ierr)
            call check_deallocate('excit_gen_pp%all_list_beta',ierr)
        end if
    
    end subroutine dealloc_excit_gen_power_pitzer_t

    subroutine dealloc_excit_gen_heat_bath_t(excit_gen_hb)

        ! Deallocate the heat bath excitation generator data.

        ! In/Out:
        !   excit_gen_hb: excit_gen_heat_bath_t to be deallocated.

        use checking, only: check_deallocate
        type(excit_gen_heat_bath_t), intent(inout) :: excit_gen_hb
        integer :: ierr

        if (allocated(excit_gen_hb%i_weights)) then
            deallocate(excit_gen_hb%i_weights, stat=ierr)
            call check_deallocate('excit_gen_hb%i_weights',ierr)
        end if
        if (allocated(excit_gen_hb%ij_weights)) then
            deallocate(excit_gen_hb%ij_weights, stat=ierr)
            call check_deallocate('excit_gen_hb%ij_weights',ierr)
        end if
        
        call dealloc_alias_table_data_t(excit_gen_hb%hb_ija)
        call dealloc_alias_table_data_t(excit_gen_hb%hb_ijab)

    end subroutine dealloc_excit_gen_heat_bath_t

!--------------------------------------------------------------------------------!
    ! alloc_alias_table_t(alias_data)
        ! Allocate alias table data.

        ! In/Out:
        !   alias_data: alias_table_data_one_ind_t to be deallocated.
    subroutine alloc_alias_table_data_one_ind_t(alias_data, len_one)
 
        use checking, only: check_allocate

        type(alias_table_data_one_ind_t), intent(inout) :: alias_data
        integer, intent(in) :: len_one
        integer :: ierr

        allocate(alias_data%weights(len_one), stat=ierr)
        call check_allocate('alias_data%weights', len_one, ierr)
        allocate(alias_data%aliasU(len_one), stat=ierr)
        call check_allocate('alias_data%aliasU', len_one, ierr)
        allocate(alias_data%aliasK(len_one), stat=ierr)
        call check_allocate('alias_data%aliasK', len_one, ierr)

    end subroutine alloc_alias_table_data_one_ind_t
    
    subroutine alloc_alias_table_data_two_ind_t(alias_data, len_two, len_one)

        use checking, only: check_allocate
        
        type(alias_table_data_two_ind_t), intent(inout) :: alias_data
        integer, intent(in) :: len_one, len_two
        integer :: ierr

        allocate(alias_data%weights(len_two, len_one), stat=ierr)
        call check_allocate('alias_data%weights', (len_one*len_two), ierr)
        allocate(alias_data%aliasU(len_two, len_one), stat=ierr)
        call check_allocate('alias_data%aliasU', (len_one*len_two), ierr)
        allocate(alias_data%aliasK(len_two, len_one), stat=ierr)
        call check_allocate('alias_data%aliasK', (len_one*len_two), ierr)
        allocate(alias_data%weights_tot(len_one), stat=ierr)
        call check_allocate('alias_data%weights_tot', len_one, ierr)

    end subroutine alloc_alias_table_data_two_ind_t

    subroutine alloc_alias_table_data_three_ind_t(alias_data, len_three, len_two, len_one)

        use checking, only: check_allocate
        
        type(alias_table_data_three_ind_t), intent(inout) :: alias_data
        integer, intent(in) :: len_one, len_two, len_three
        integer :: ierr

        allocate(alias_data%weights(len_three, len_two, len_one), stat=ierr)
        call check_allocate('alias_data%weights', (len_one*len_two*len_three), ierr)
        allocate(alias_data%aliasU(len_three, len_two, len_one), stat=ierr)
        call check_allocate('alias_data%aliasU', (len_one*len_two*len_three), ierr)
        allocate(alias_data%aliasK(len_three, len_two, len_one), stat=ierr)
        call check_allocate('alias_data%aliasK', (len_one*len_two*len_three), ierr)
        allocate(alias_data%weights_tot(len_two, len_one), stat=ierr)
        call check_allocate('alias_data%weights_tot', (len_one*len_two), ierr)

    end subroutine alloc_alias_table_data_three_ind_t
    
    subroutine alloc_alias_table_data_2indarray_three_ind_t(alias_data, len_three, array_two, len_one)

        use checking, only: check_allocate
        
        type(alias_table_data_three_ind_t), intent(inout) :: alias_data
        integer, intent(in) :: len_one, array_two(2), len_three
        integer :: ierr, len_two

        len_two = array_two(2) - array_two(1) + 1

        allocate(alias_data%weights(len_three, array_two(1):array_two(2), len_one), stat=ierr)
        call check_allocate('alias_data%weights', (len_one*len_two*len_three), ierr)
        allocate(alias_data%aliasU(len_three, array_two(1):array_two(2), len_one), stat=ierr)
        call check_allocate('alias_data%aliasU', (len_one*len_two*len_three), ierr)
        allocate(alias_data%aliasK(len_three, array_two(1):array_two(2), len_one), stat=ierr)
        call check_allocate('alias_data%aliasK', (len_one*len_two*len_three), ierr)
        allocate(alias_data%weights_tot(array_two(1):array_two(2), len_one), stat=ierr)
        call check_allocate('alias_data%weights_tot', (len_one*len_two), ierr)

    end subroutine alloc_alias_table_data_2indarray_three_ind_t

    subroutine alloc_alias_table_data_four_ind_t(alias_data, len_four, len_three, len_two, len_one)

        use checking, only: check_allocate
        
        type(alias_table_data_four_ind_t), intent(inout) :: alias_data
        integer, intent(in) :: len_one, len_two, len_three, len_four
        integer :: ierr

        allocate(alias_data%weights(len_four, len_three, len_two, len_one), stat=ierr)
        call check_allocate('alias_data%weights', (len_one*len_two*len_three*len_four), ierr)
        allocate(alias_data%aliasU(len_four, len_three, len_two, len_one), stat=ierr)
        call check_allocate('alias_data%aliasU', (len_one*len_two*len_three*len_four), ierr)
        allocate(alias_data%aliasK(len_four, len_three, len_two, len_one), stat=ierr)
        call check_allocate('alias_data%aliasK', (len_one*len_two*len_three*len_four), ierr)
        allocate(alias_data%weights_tot(len_three, len_two, len_one), stat=ierr)
        call check_allocate('alias_data%weights_tot', (len_one*len_two*len_three), ierr)

    end subroutine alloc_alias_table_data_four_ind_t
!--------------------------------------------------------------------------------!

!--------------------------------------------------------------------------------!
    ! dealloc_alias_table_t(alias_data)
        ! Deallocate alias table data.

        ! In/Out:
        !   alias_data: alias_table_data_one_ind_t to be deallocated.
    subroutine dealloc_alias_table_data_one_ind_t(alias_data)
        
        use checking, only: check_deallocate
        type(alias_table_data_one_ind_t), intent(inout) :: alias_data
        integer :: ierr

        if (allocated(alias_data%weights)) then
            deallocate(alias_data%weights, stat=ierr)
            call check_deallocate('alias_data%weights',ierr)
        end if
        if (allocated(alias_data%aliasU)) then
            deallocate(alias_data%aliasU, stat=ierr)
            call check_deallocate('alias_data%aliasU',ierr)
        end if
        if (allocated(alias_data%aliasK)) then
            deallocate(alias_data%aliasK, stat=ierr)
            call check_deallocate('alias_data%aliasK',ierr)
        end if

    end subroutine dealloc_alias_table_data_one_ind_t
    
    subroutine dealloc_alias_table_data_two_ind_t(alias_data)

        use checking, only: check_deallocate
        type(alias_table_data_two_ind_t), intent(inout) :: alias_data
        integer :: ierr

        if (allocated(alias_data%weights)) then
            deallocate(alias_data%weights, stat=ierr)
            call check_deallocate('alias_data%weights',ierr)
        end if
        if (allocated(alias_data%aliasU)) then
            deallocate(alias_data%aliasU, stat=ierr)
            call check_deallocate('alias_data%aliasU',ierr)
        end if
        if (allocated(alias_data%aliasK)) then
            deallocate(alias_data%aliasK, stat=ierr)
            call check_deallocate('alias_data%aliasK',ierr)
        end if
        if (allocated(alias_data%weights_tot)) then
            deallocate(alias_data%weights_tot, stat=ierr)
            call check_deallocate('alias_data%weights_tot',ierr)
        end if

    end subroutine dealloc_alias_table_data_two_ind_t

    subroutine dealloc_alias_table_data_three_ind_t(alias_data)

        use checking, only: check_deallocate
        type(alias_table_data_three_ind_t), intent(inout) :: alias_data
        integer :: ierr

        if (allocated(alias_data%weights)) then
            deallocate(alias_data%weights, stat=ierr)
            call check_deallocate('alias_data%weights',ierr)
        end if
        if (allocated(alias_data%aliasU)) then
            deallocate(alias_data%aliasU, stat=ierr)
            call check_deallocate('alias_data%aliasU',ierr)
        end if
        if (allocated(alias_data%aliasK)) then
            deallocate(alias_data%aliasK, stat=ierr)
            call check_deallocate('alias_data%aliasK',ierr)
        end if
        if (allocated(alias_data%weights_tot)) then
            deallocate(alias_data%weights_tot, stat=ierr)
            call check_deallocate('alias_data%weights_tot',ierr)
        end if

    end subroutine dealloc_alias_table_data_three_ind_t

    subroutine dealloc_alias_table_data_four_ind_t(alias_data)

        use checking, only: check_deallocate
        type(alias_table_data_four_ind_t), intent(inout) :: alias_data
        integer :: ierr

        if (allocated(alias_data%weights)) then
            deallocate(alias_data%weights, stat=ierr)
            call check_deallocate('alias_data%weights',ierr)
        end if
        if (allocated(alias_data%aliasU)) then
            deallocate(alias_data%aliasU, stat=ierr)
            call check_deallocate('alias_data%aliasU',ierr)
        end if
        if (allocated(alias_data%aliasK)) then
            deallocate(alias_data%aliasK, stat=ierr)
            call check_deallocate('alias_data%aliasK',ierr)
        end if
        if (allocated(alias_data%weights_tot)) then
            deallocate(alias_data%weights_tot, stat=ierr)
            call check_deallocate('alias_data%weights_tot',ierr)
        end if

    end subroutine dealloc_alias_table_data_four_ind_t
!--------------------------------------------------------------------------------!
end module excit_gens
