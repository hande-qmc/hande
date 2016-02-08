module excit_gens
use const
implicit none

!Type containing data for excitation generators
type excit_gen_data_t
    ! Probability of attempting single or double excitations.
    real(p) :: pattempt_single, pattempt_double
end type excit_gen_data_t

end module excit_gens
