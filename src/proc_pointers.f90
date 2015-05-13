module proc_pointers

use const, only: i0, p, dp, int_64, int_p
use determinants, only: det_info_t
use excitations, only: excit_t

implicit none

! Note that Intel compilers (at least 11.1) don't like having using a variable
! that's imported as a array size in abstract interfaces.

abstract interface
    pure subroutine i_decoder(sys,f,d)
        use system, only: sys_t
        import :: i0, det_info_t
        implicit none
        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%string_len)
        type(det_info_t), intent(inout) :: d
    end subroutine i_decoder
    pure subroutine i_update_proj_energy(sys, f0, d, pop, D0_pop_sum, proj_energy_sum, excitation, hmatel)
        use system, only: sys_t
        import :: det_info_t, p, i0, excit_t
        implicit none
        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f0(:)
        type(det_info_t), intent(in) :: d
        real(p), intent(in) :: pop
        real(p), intent(inout) :: D0_pop_sum, proj_energy_sum
        type(excit_t), intent(inout) :: excitation
        real(p), intent(out) :: hmatel
    end subroutine i_update_proj_energy
    subroutine i_update_proj_hfs(sys, f, fpop, f_hfpop, fdata, excitation, hmatel,&
                                     D0_hf_pop,proj_hf_O_hpsip, proj_hf_H_hfpsip)
        use system, only: sys_t
        import :: i0, p, excit_t
        implicit none
        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(:)
        integer, intent(in) :: fpop, f_hfpop
        real(p), intent(in) :: fdata(:), hmatel
        type(excit_t), intent(in) :: excitation
        real(p), intent(inout) :: D0_hf_pop, proj_hf_O_hpsip, proj_hf_H_hfpsip
    end subroutine i_update_proj_hfs
    subroutine i_update_dmqmc_energy_and_trace(sys, excitation, d, H00, walker_pop, diag, trace, energy)
        use system, only: sys_t
        import :: excit_t, p, det_info_t
        implicit none
        type(sys_t), intent(in) :: sys
        type(excit_t), intent(inout) :: excitation
        type(det_info_t), intent(in) :: d
        real(p), intent(in) :: H00, walker_pop
        real(p), intent(in) :: diag
        real(p), intent(inout) :: trace(:)
        real(p), intent(inout) :: energy
    end subroutine i_update_dmqmc_energy_and_trace
    subroutine i_update_dmqmc_estimators(sys, cdet, excitation, H00, walker_pop, estimate)
        use system, only: sys_t
        use determinants, only: det_info_t
        import :: excit_t, p
        implicit none
        type(sys_t), intent(in) :: sys
        type(det_info_t), intent(in) :: cdet
        type(excit_t), intent(in) :: excitation
        real(p), intent(in) :: H00, walker_pop
        real(p), intent(inout) :: estimate
    end subroutine i_update_dmqmc_estimators
    subroutine i_update_dmqmc_correlation_function(sys, cdet, excitation, H00, walker_pop, mask, cfunc)
        use system, only: sys_t
        use determinants, only: det_info_t
        import :: excit_t, p, i0
        implicit none
        type(sys_t), intent(in) :: sys
        type(det_info_t), intent(in) :: cdet
        type(excit_t), intent(in) :: excitation
        real(p), intent(in) :: H00, walker_pop
        integer(i0), allocatable, intent(in) :: mask(:)
        real(p), intent(inout) :: cfunc
    end subroutine i_update_dmqmc_correlation_function
    subroutine i_gen_excit(rng, sys, qmc_in, d, pgen, connection, hmatel)
        use dSFMT_interface, only: dSFMT_t
        use qmc_data, only: qmc_in_t
        use system, only: sys_t
        import :: det_info_t, excit_t, p
        implicit none
        type(dSFMT_t), intent(inout) :: rng
        type(sys_t), intent(in) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        type(det_info_t), intent(in) :: d
        real(p), intent(out) :: pgen, hmatel
        type(excit_t), intent(out) :: connection
    end subroutine i_gen_excit
    subroutine i_gen_excit_finalise(rng, sys, qmc_in, d, connection, hmatel)
        use dSFMT_interface, only: dSFMT_t
        use qmc_data, only: qmc_in_t
        use system, only: sys_t
        import :: det_info_t, excit_t, p
        implicit none
        type(dSFMT_t), intent(inout) :: rng
        type(sys_t), intent(in) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        type(det_info_t), intent(in) :: d
        type(excit_t), intent(inout) :: connection
        real(p), intent(out) :: hmatel
    end subroutine i_gen_excit_finalise
    pure function i_sc0(sys, f) result(hmatel)
        use system, only: sys_t
        import :: p, i0
        implicit none
        real(p) :: hmatel
        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%string_len)
    end function i_sc0
    pure function i_ex0(sys, occ_list) result(hmatel)
        use system, only: sys_t
        import :: p
        implicit none
        real(p) :: hmatel
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: occ_list(:)
    end function i_ex0
    subroutine i_create_spawned_particle(basis, reference, d, connection, nspawned, particle_indx, spawn, f)
        use basis_types, only: basis_t
        use spawn_data, only: spawn_t
        use qmc_data, only: reference_t
        import :: excit_t, det_info_t, int_p, i0
        implicit none
        type(basis_t), intent(in) :: basis
        type(reference_t), intent(in) :: reference
        type(det_info_t), intent(in) :: d
        type(excit_t), intent(in) :: connection
        integer(int_p), intent(in) :: nspawned
        integer, intent(in) :: particle_indx
        integer(i0), intent(in), optional, target :: f(:)
        type(spawn_t), intent(inout) :: spawn
    end subroutine i_create_spawned_particle
    subroutine i_create_spawned_particle_dm(basis, reference, f1, f2, connection, nspawned, spawning_end, particle_indx, spawn)
        use spawn_data, only: spawn_t
        use basis_types, only: basis_t
        use qmc_data, only: reference_t
        import :: excit_t, i0, int_p
        implicit none
        type(basis_t), intent(in) :: basis
        type(reference_t), intent(in) :: reference
        integer(i0), intent(in) :: f1(basis%string_len)
        integer(i0), intent(in) :: f2(basis%string_len)
        type(excit_t), intent(in) :: connection
        integer(int_p), intent(in) :: nspawned
        integer, intent(in) :: spawning_end, particle_indx
        type(spawn_t), intent(inout) :: spawn
    end subroutine i_create_spawned_particle_dm
    subroutine i_trial_fn(sys, cdet, connection, weights, hmatel)
        use system, only: sys_t
        import :: det_info_t, excit_t, p
        type(sys_t), intent(in) :: sys
        type(det_info_t), intent(in) :: cdet
        type(excit_t), intent(in) :: connection
        real(p), allocatable, intent(in) :: weights(:)
        real(p), intent(inout) :: hmatel
    end subroutine i_trial_fn

    ! generic procedures...
    subroutine i_sub()
    end subroutine i_sub

end interface

procedure(i_decoder), pointer :: decoder_ptr => null()

procedure(i_update_proj_energy), pointer :: update_proj_energy_ptr => null()
procedure(i_update_proj_hfs), pointer :: update_proj_hfs_ptr => null()

procedure(i_update_dmqmc_energy_and_trace), pointer :: update_dmqmc_energy_and_trace_ptr => null()
procedure(i_update_dmqmc_estimators), pointer :: update_dmqmc_energy_squared_ptr => null()
procedure(i_update_dmqmc_estimators), pointer :: update_dmqmc_stag_mag_ptr => null()
procedure(i_update_dmqmc_correlation_function), pointer :: update_dmqmc_correlation_ptr => null()

procedure(i_sc0), pointer :: sc0_ptr => null()
procedure(i_sc0), pointer :: op0_ptr => null()
procedure(i_sc0), pointer :: trial_dm_ptr => null()
! [review] - JSS: is this really on a performance-critical path to be worth the hassle of making it a procedure pointer?
procedure(i_ex0), pointer :: ex0_ptr => null()

procedure(i_sub), pointer :: dmqmc_initial_distribution_ptr => null()

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
    subroutine i_spawner(rng, sys, qmc_in, tau, spawn_cutoff, real_factor, d, parent_sign, gen_excit_ptr, weights, &
                         nspawned, connection)
        use dSFMT_interface, only: dSFMT_t
        use qmc_data, only: qmc_in_t
        use system, only: sys_t
        import :: det_info_t, excit_t, gen_excit_ptr_t, int_p, p, dp
        implicit none
        type(dSFMT_t), intent(inout) :: rng
        type(sys_t), intent(in) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        real(p), intent(in) :: tau
        integer(int_p), intent(in) :: spawn_cutoff
        integer(int_p), intent(in) :: real_factor
        type(det_info_t), intent(in) :: d
        integer(int_p), intent(in) :: parent_sign
        type(gen_excit_ptr_t), intent(in) :: gen_excit_ptr
        real(p), allocatable, intent(in) :: weights(:)
        integer(int_p), intent(out) :: nspawned
        type(excit_t), intent(out) :: connection
    end subroutine i_spawner
end interface

procedure(i_spawner), pointer :: spawner_ptr => null()
procedure(i_spawner), pointer :: spawner_hfs_ptr => null()

end module proc_pointers
