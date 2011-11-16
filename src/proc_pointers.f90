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
    subroutine i_update_dmqmc_energy(idet,beta_index)
        implicit none
        integer, intent(in) :: idet, beta_index
    end subroutine i_update_dmqmc_energy
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
    subroutine i_set_parent_flag(pop, f, flag)
        use basis, only: basis_length
        import :: i0
        implicit none
        integer, intent(in) :: pop
        integer(i0), intent(in) :: f(basis_length)
        integer, intent(out) :: flag
    end subroutine i_set_parent_flag
    subroutine i_create_spawned_particle(d, connection, nspawned, spawned_pop)
        import :: excit, det_info 
        implicit none
        type(det_info), intent(in) :: d
        type(excit), intent(in) :: connection
        integer, intent(in) :: nspawned, spawned_pop
    end subroutine i_create_spawned_particle
    subroutine i_create_spawned_particle_dm(f1, f2, connection, nspawned, spawning_end)
        use basis, only: basis_length
        import :: excit, i0
        implicit none
        integer(i0), intent(in) :: f1(basis_length)
        integer(i0), intent(in) :: f2(basis_length)
        type(excit), intent(in) :: connection
        integer, intent(in) :: nspawned, spawning_end
    end subroutine i_create_spawned_particle_dm

    ! generic procedures...
    subroutine i_sub()
    end subroutine i_sub

end interface

procedure(i_decoder), pointer :: decoder_ptr => null()
procedure(i_update_proj_energy), pointer :: update_proj_energy_ptr => null()
procedure(i_update_dmqmc_energy), pointer :: update_dmqmc_energy_ptr => null()
procedure(i_spawner), pointer :: spawner_ptr => null()
procedure(i_sc0), pointer :: sc0_ptr => null()
procedure(i_sub), pointer :: annihilate_main_list_ptr => null()
procedure(i_sub), pointer :: annihilate_spawned_list_ptr => null()
procedure(i_set_parent_flag), pointer :: set_parent_flag_ptr => null()
procedure(i_create_spawned_particle), pointer :: create_spawned_particle_ptr => null()
procedure(i_create_spawned_particle_dm), pointer :: create_spawned_particle_dm_ptr => null()

end module proc_pointers
