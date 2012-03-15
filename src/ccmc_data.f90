module ccmc_data

use const, only: i0, p

implicit none

! Work around Fortran's lack of arrays of pointers...
type bit_string_ptr
    integer(i0), pointer :: f(:)
end type bit_string_ptr

! Information about a cluster of excitors in addition to that stored in
! a det_info variable.
type cluster_t
    ! Pointers to the determinants formed by applying individual excitors to the
    ! reference determinant.
    type(bit_string_ptr), allocatable :: excitors(:) ! max: truncation_level+2
    ! Number of excitors in cluster
    integer :: nexcitors
    ! Excitation level relative to the reference determinant of the determinant
    ! formed by applying the cluster of excitors to the reference determinant.
    integer :: excitation_level
    ! Overall amplitude of the cluster.
    real(p) :: amplitude
    ! < D_i | a_i D_0 >, where D_i is the determinant formed by applying the
    ! cluster of excitors to the reference determinant.  Equal to +1 or -1.
    integer :: cluster_to_det_sign
    ! probability of selecting this cluster at random
    real(p) :: pselect
end type cluster_t

end module ccmc_data
