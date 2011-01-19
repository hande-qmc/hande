module proc_pointers

use const, only: i0, p
use determinants, only: det_info
use excitations, only: excit

implicit none

! Note that Intel compilers (at least 11.1) don't like having using a variable
! that's imported as a array size in abstract interfaces.

abstract interface
    subroutine i_decoder(f,d)
        use basis, only: basis_length
        import :: i0, det_info
        implicit none
        integer(i0), intent(in) :: f(basis_length)
        type(det_info), intent(inout) :: d
    end subroutine i_decoder
    subroutine i_update_proj_energy(idet)
        implicit none
        integer, intent(in) :: idet
    end subroutine i_update_proj_energy
    subroutine i_spawner(d, parent_sign, nspawned, connection)
        import :: det_info, excit
        implicit none
        type(det_info), intent(in) :: d
        integer, intent(in) :: parent_sign
        integer, intent(out) :: nspawned
        type(excit), intent(out) :: connection
    end subroutine i_spawner
    function i_sc0(f) result(hmatel)
        use basis, only: basis_length
        import :: p, i0
        implicit none
        real(p) :: hmatel
        integer(i0), intent(in) :: f(basis_length)
    end function i_sc0
end interface

procedure(i_decoder), pointer :: decoder => null()
procedure(i_update_proj_energy), pointer :: update_proj_energy => null()
procedure(i_spawner), pointer :: spawner => null()
procedure(i_sc0), pointer :: sc0 => null()

end module proc_pointers
