module excit_gens
use const
implicit none

!Type containing data for excitation generators
type excit_gen_data_t
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
    ! (0:string_len,-N:N,0,0) (1D)
    ! (0:string_len,-N:N,-N:N,0) (2D)
    ! (0:string_len,-N:N,-N:N,-N:N) (3D)
    ! NOTE: this contains values of k_i+k_j which cannot be formed by the basis with
    ! the energy cutoff.  Memory can be saved by not using a cubic array for
    ! k_i+k_j...
    integer(i0), allocatable :: ueg_ternary_conserve(:,:,:,:)
end type excit_gen_data_t

end module excit_gens
