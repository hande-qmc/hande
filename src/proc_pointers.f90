module proc_pointers

use const, only: i0, p, dp
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
    subroutine i_gen_excit(d, pgen, connection,hmatel)
        import :: det_info, excit, p
        implicit none
        type(det_info), intent(in) :: d
        real(p), intent(out) :: pgen, hmatel
        type(excit), intent(out) :: connection
    end subroutine i_gen_excit
    subroutine i_death(mat, pop, tot_pop, ndeath)
        import :: p
        implicit none
        real(p), intent(in) :: mat
        integer, intent(inout) :: pop, tot_pop
        integer, intent(out) :: ndeath
    end subroutine i_death
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

    ! generic procedures...
    subroutine i_sub()
    end subroutine i_sub
 
    !fsfciqmc related procedures......................................................................
    
    function i_rng() result(r)
        import:: dp
        implicit none
        real(dp) :: r
    end function i_rng

    !...............................................................................................

end interface

procedure(i_decoder), pointer :: decoder_ptr => null()
procedure(i_update_proj_energy), pointer :: update_proj_energy_ptr => null()
procedure(i_spawner), pointer :: spawner_ptr => null()
procedure(i_gen_excit), pointer :: gen_excit_ptr => null()
procedure(i_death), pointer :: death_ptr => null()
procedure(i_sc0), pointer :: sc0_ptr => null()
procedure(i_sub), pointer :: annihilate_main_list_ptr => null()
procedure(i_sub), pointer :: annihilate_spawned_list_ptr => null()
procedure(i_set_parent_flag), pointer :: set_parent_flag_ptr => null()
procedure(i_create_spawned_particle), pointer :: create_spawned_particle_ptr => null()

!fsfciqmc related procedures......................................................................
procedure(i_rng), pointer :: rng_ptr => null() 
procedure(i_sc0), pointer :: system_sc0_ptr => null()
procedure(i_decoder), pointer :: system_decoder_ptr => null()
procedure(i_update_proj_energy), pointer :: system_update_proj_energy_ptr => null()
    


!...............................................................................................

end module proc_pointers
