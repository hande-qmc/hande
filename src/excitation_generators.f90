module excitations

! Module for generating random excitations.

use const

implicit none

contains

    function excit_gen_norm_hub_k(f) result(norm)

        use basis, only: basis_length

        real(dp) :: norm
        integer(i0), intent(in) :: f(basis_length)

    end function excit_gen_norm_hub_k

end module excitations
