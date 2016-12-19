module excit_gens
use const
implicit none

!Data for the Cauchy_Schwarz excit gen

! [review] - JSS: unclear why a new named parameter for int_32.
!The integer types have been chosen to be int_32
integer(int_32), parameter :: int_bas = int_32

! [review] - JSS: document what each component holds.
type excit_gen_cauchy_schwarz_t
    real(p), allocatable :: ia_aliasP(:,:) !(max(sys%nvirt_alpha,sys%nvirt_beta),sys%nel)
    integer(int_bas), allocatable :: ia_aliasY(:,:) !(max(sys%nvirt_alpha,sys%nvirt_beta),sys%nel)
    real(p), allocatable :: ia_weights(:,:) !(max(sys%nvirt_alpha,sys%nvirt_beta),sys%nel)
    real(p), allocatable :: ia_weights_tot(:) !(sys%nel)
    real(p), allocatable :: jb_aliasP(:,:,:) !(maxval(sys%read_in%pg_sym%nbasis_sym_spin),sys%sym0_tot:sys%sym_max_tot,sys%nel)
    integer(int_bas), allocatable :: jb_aliasY(:,:,:) !(maxval(sys%read_in%pg_sym%nbasis_sym_spin),sys%sym0_tot:sys%sym_max_tot,sys%nel)
    real(p), allocatable :: jb_weights(:,:,:) !(maxval(sys%read_in%pg_sym%nbasis_sym_spin),sys%sym0_tot:sys%sym_max_tot,sys%nel)
    real(p), allocatable :: jb_weights_tot(:,:) !(sys%sym0_tot:sys%sym_max_tot,sys%nel)
    integer(int_bas), allocatable :: virt_list_alpha(:) !(sys%nvirt_alpha)
    integer(int_bas), allocatable :: virt_list_beta(:) !(sys%nvirt_beta)
    integer(int_bas), allocatable :: occ_list(:) !(sys%nel+1)  !The +1 is a pad
end type excit_gen_cauchy_schwarz_t


!Type containing data for excitation generators
type excit_gen_data_t
    ! Excitation generator to use, duplicated from qmc_in.
    integer :: excit_gen

    ! Probability of attempting single or double excitations.
    real(p) :: pattempt_single, pattempt_double

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
    type(excit_gen_cauchy_schwarz_t) :: excit_gen_cs
end type excit_gen_data_t

contains

    subroutine dealloc_excit_gen_data_t(excit_gen_data)

        ! Deallocate the excitation generator data.

        ! In/Out:
        !   excit_gen_data: excit_gen_data_t to be deallocated.

        ! It is a bit superfluous to have this subroutine for one deallocation, but it will be more
        ! useful in future.

        type(excit_gen_data_t), intent(inout) :: excit_gen_data

        if (allocated(excit_gen_data%ueg_ternary_conserve)) deallocate(excit_gen_data%ueg_ternary_conserve)

        call dealloc_excit_gen_cauchy_schwarz_t(excit_gen_data%excit_gen_cs)

    end subroutine dealloc_excit_gen_data_t

    subroutine dealloc_excit_gen_cauchy_schwarz_t(excit_gen_cs)
        
        ! Deallocate the cauchy schwarz excitation generator data.

        ! In/Out:
        !   excit_gen_cs: excit_gen_cauchy_schwarz_t to be deallocated.

        type(excit_gen_cauchy_schwarz_t), intent(inout) :: excit_gen_cs

        if (allocated(excit_gen_cs%ia_aliasP)) deallocate(excit_gen_cs%ia_aliasP)
        if (allocated(excit_gen_cs%ia_aliasY)) deallocate(excit_gen_cs%ia_aliasY)
        if (allocated(excit_gen_cs%ia_weights)) deallocate(excit_gen_cs%ia_weights)
        if (allocated(excit_gen_cs%ia_weights_tot)) deallocate(excit_gen_cs%ia_weights_tot)
        if (allocated(excit_gen_cs%jb_aliasP)) deallocate(excit_gen_cs%jb_aliasP)
        if (allocated(excit_gen_cs%jb_aliasY)) deallocate(excit_gen_cs%jb_aliasY)
        if (allocated(excit_gen_cs%jb_weights)) deallocate(excit_gen_cs%jb_weights)
        if (allocated(excit_gen_cs%jb_weights_tot)) deallocate(excit_gen_cs%jb_weights_tot)
        if (allocated(excit_gen_cs%virt_list_alpha)) deallocate(excit_gen_cs%virt_list_alpha)
        if (allocated(excit_gen_cs%virt_list_beta)) deallocate(excit_gen_cs%virt_list_beta)
        if (allocated(excit_gen_cs%occ_list)) deallocate(excit_gen_cs%occ_list)

    end subroutine dealloc_excit_gen_cauchy_schwarz_t

end module excit_gens
