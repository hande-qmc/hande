module proc_pointers

use const, only: i0, p, dp, lint
use determinants, only: det_info
use excitations, only: excit

implicit none

! Note that Intel compilers (at least 11.1) don't like having using a variable
! that's imported as a array size in abstract interfaces.

abstract interface
    pure subroutine i_decoder(f,d)
        use basis, only: basis_length
        import :: i0, det_info
        implicit none
        integer(i0), intent(in) :: f(basis_length)
        type(det_info), intent(inout) :: d
    end subroutine i_decoder
    pure subroutine i_update_proj_energy(f0, d, pop, D0_pop_sum, proj_energy_sum, excitation, hmatel)
        import :: det_info, p, i0, excit
        implicit none
        integer(i0), intent(in) :: f0(:)
        type(det_info), intent(in) :: d
        real(p), intent(in) :: pop
        real(p), intent(inout) :: D0_pop_sum, proj_energy_sum
        type(excit), intent(out) :: excitation
        real(p), intent(out) :: hmatel
    end subroutine i_update_proj_energy
    subroutine i_update_proj_hfs(f, fpop, f_hfpop, fdata, excitation, hmatel,&
                                     D0_hf_pop,proj_hf_O_hpsip, proj_hf_H_hfpsip)
        use basis, only: basis_length
        import :: i0, p, excit
        implicit none
        integer(i0), intent(in) :: f(basis_length)
        integer, intent(in) :: fpop, f_hfpop
        real(p), intent(in) :: fdata(:), hmatel
        type(excit), intent(in) :: excitation
        real(p), intent(inout) :: D0_hf_pop, proj_hf_O_hpsip, proj_hf_H_hfpsip
    end subroutine i_update_proj_hfs
    subroutine i_update_dmqmc_estimators(idet,excitation,walker_pop)
        import :: excit, p
        implicit none
        integer, intent(in) :: idet
        type(excit), intent(in) :: excitation
        real(p), intent(in) :: walker_pop
    end subroutine i_update_dmqmc_estimators
    subroutine i_gen_excit(rng, d, pgen, connection, hmatel)
        use dSFMT_interface, only: dSFMT_t
        import :: det_info, excit, p
        implicit none
        type(dSFMT_t), intent(inout) :: rng
        type(det_info), intent(in) :: d
        real(p), intent(out) :: pgen, hmatel
        type(excit), intent(out) :: connection
    end subroutine i_gen_excit
    subroutine i_gen_excit_finalise(rng, d, connection, hmatel)
        use dSFMT_interface, only: dSFMT_t
        import :: det_info, excit, p
        implicit none
        type(dSFMT_t), intent(inout) :: rng
        type(det_info), intent(in) :: d
        type(excit), intent(inout) :: connection
        real(p), intent(out) :: hmatel
    end subroutine i_gen_excit_finalise
    subroutine i_death(rng, mat, pop, tot_pop, ndeath)
        use dSFMT_interface, only: dSFMT_t
        import :: p, lint
        implicit none
        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(in) :: mat
        integer, intent(inout) :: pop, ndeath
        integer(lint), intent(inout) :: tot_pop
    end subroutine i_death
    pure function i_sc0(f) result(hmatel)
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
    subroutine i_create_spawned_particle(d, connection, nspawned, particle_indx, spawn)
        use spawn_data, only: spawn_t
        import :: excit, det_info
        implicit none
        type(det_info), intent(in) :: d
        type(excit), intent(in) :: connection
        integer, intent(in) :: nspawned, particle_indx
        type(spawn_t), intent(inout) :: spawn
    end subroutine i_create_spawned_particle
    subroutine i_create_spawned_particle_dm(f1, f2, connection, nspawned, spawning_end, particle_indx, spawn)
        use basis, only: basis_length
        use spawn_data, only: spawn_t
        import :: excit, i0
        implicit none
        integer(i0), intent(in) :: f1(basis_length)
        integer(i0), intent(in) :: f2(basis_length)
        type(excit), intent(in) :: connection
        integer, intent(in) :: nspawned, spawning_end, particle_indx
        type(spawn_t), intent(inout) :: spawn
    end subroutine i_create_spawned_particle_dm
    subroutine i_trial_fn(cdet, connection, hmatel)
        import :: det_info, excit, p
        type(det_info), intent(in) :: cdet
        type(excit), intent(in) :: connection
        real(p), intent(inout) :: hmatel
    end subroutine i_trial_fn

    ! generic procedures...
    subroutine i_sub()
    end subroutine i_sub

end interface

procedure(i_decoder), pointer :: decoder_ptr => null()

procedure(i_update_proj_energy), pointer :: update_proj_energy_ptr => null()
procedure(i_update_proj_hfs), pointer :: update_proj_hfs_ptr => null()

procedure(i_update_dmqmc_estimators), pointer :: update_dmqmc_energy_ptr => null()
procedure(i_update_dmqmc_estimators), pointer :: update_dmqmc_energy_squared_ptr => null()
procedure(i_update_dmqmc_estimators), pointer :: update_dmqmc_stag_mag_ptr => null()
procedure(i_update_dmqmc_estimators), pointer :: update_dmqmc_correlation_ptr => null()

procedure(i_death), pointer :: death_ptr => null()

procedure(i_sc0), pointer :: sc0_ptr => null()
procedure(i_sc0), pointer :: op0_ptr => null()

procedure(i_sub), pointer :: dmqmc_initial_distribution_ptr => null()

procedure(i_set_parent_flag), pointer :: set_parent_flag_ptr => null()

procedure(i_create_spawned_particle), pointer :: create_spawned_particle_ptr => null()
procedure(i_create_spawned_particle_dm), pointer :: create_spawned_particle_dm_ptr => null()


! Single structure for all types of excitation generator so we can use the same
! interface for spawning routines which use different types of generator.
type gen_excit_ptr_t
    procedure(i_gen_excit), nopass, pointer :: full => null()
    procedure(i_gen_excit), nopass, pointer :: init => null()
    procedure(i_gen_excit_finalise), nopass, pointer :: finalise => null()
    procedure(i_trial_fn), nopass, pointer :: trial_fn => null()
end type gen_excit_ptr_t

type(gen_excit_ptr_t) :: gen_excit_ptr, gen_excit_hfs_ptr

abstract interface
    subroutine i_spawner(rng, d, parent_sign, gen_excit_ptr, nspawned, connection)
        use dSFMT_interface, only: dSFMT_t
        import :: det_info, excit, gen_excit_ptr_t
        implicit none
        type(dSFMT_t), intent(inout) :: rng
        type(det_info), intent(in) :: d
        integer, intent(in) :: parent_sign
        type(gen_excit_ptr_t), intent(in) :: gen_excit_ptr
        integer, intent(out) :: nspawned
        type(excit), intent(out) :: connection
    end subroutine i_spawner
end interface

procedure(i_spawner), pointer :: spawner_ptr => null()
procedure(i_spawner), pointer :: spawner_hfs_ptr => null()

end module proc_pointers
