module dealloc

! Functions to deallocate various derived types.

implicit none

contains

    subroutine dealloc_qmc_state_t(qs)

        ! Deallocate all components of a qmc_state_t object

        ! In/Out:
        !   qs: qmc_state to be deallocated.

        use qmc_data, only: qmc_state_t
        use particle_t_utils, only: dealloc_particle_t
        use spawn_data, only: dealloc_spawn_t
        use reference_determinant, only: dealloc_reference_t
        use load_balancing, only: dealloc_parallel_t

        type(qmc_state_t), intent(inout) :: qs

        if (allocated(qs%shift)) deallocate(qs%shift)
        if (allocated(qs%vary_shift)) deallocate(qs%vary_shift)

        if (allocated(qs%trial%wfn_dat)) deallocate(qs%trial%wfn_dat)

        if (allocated(qs%excit_gen_data%ueg_ternary_conserve)) deallocate(qs%excit_gen_data%ueg_ternary_conserve)

        call dealloc_particle_t(qs%psip_list)
        call dealloc_spawn_t(qs%spawn_store%spawn)
        call dealloc_spawn_t(qs%spawn_store%spawn_recv)
        call dealloc_reference_t(qs%ref)
        call dealloc_parallel_t(qs%par_info)

    end subroutine dealloc_qmc_state_t

    subroutine dealloc_sys_t(sys)

        ! Deallocate all allocated components of a sys_t object

        ! In/Out:
        !   sys: sys_t object to be deallocated.

        use system, only: sys_t
        use basis_types, only: dealloc_basis_t

        type(sys_t), intent(inout) :: sys

        call dealloc_basis_t(sys%basis)
        call dealloc_sys_lattice_t(sys%lattice)
        call dealloc_sys_k_lattice_t(sys%k_lattice)
        call dealloc_sys_real_lattice_t(sys%real_lattice)
        call dealloc_sys_hubbard_t(sys%hubbard)
        call dealloc_sys_heisenberg_t(sys%heisenberg)
        call dealloc_sys_ueg_t(sys%ueg)
        call dealloc_sys_read_in_t(sys%read_in)

    end subroutine dealloc_sys_t

    subroutine dealloc_sys_lattice_t(lattice)

        ! Deallocate all allocated components of a sys_lattice_t object

        ! In/Out:
        !   sys_lattice: sys_lattice_t object to be deallocated.

        use system, only: sys_lattice_t

        type(sys_lattice_t), intent(inout) :: lattice

        if (allocated(lattice%lattice)) deallocate(lattice%lattice)
        if (allocated(lattice%rlattice)) deallocate(lattice%rlattice)
        if (allocated(lattice%box_length)) deallocate(lattice%box_length)

    end subroutine dealloc_sys_lattice_t

    subroutine dealloc_sys_k_lattice_t(k_lattice)

        ! Deallocate all allocated components of a sys_k_lattice_t object

        ! In/Out:
        !   sys_k_lattice: sys_k_lattice_t object to be deallocated.

        use system, only: sys_k_lattice_t

        type(sys_k_lattice_t), intent(inout) :: k_lattice

        if (allocated(k_lattice%ktwist)) deallocate(k_lattice%ktwist)

    end subroutine dealloc_sys_k_lattice_t

    subroutine dealloc_sys_real_lattice_t(real_lattice)

        ! Deallocate all allocated components of a sys_real_lattice_t object

        ! In/Out:
        !   sys_real_lattice: sys_real_lattice_t object to be deallocated.

        use system, only: sys_real_lattice_t

        type(sys_real_lattice_t), intent(inout) :: real_lattice

        if (allocated(real_lattice%tmat)) deallocate(real_lattice%tmat)
        if (allocated(real_lattice%connected_orbs)) deallocate(real_lattice%connected_orbs)
        if (allocated(real_lattice%connected_sites)) deallocate(real_lattice%connected_sites)
        if (allocated(real_lattice%next_nearest_orbs)) deallocate(real_lattice%next_nearest_orbs)

    end subroutine dealloc_sys_real_lattice_t

    subroutine dealloc_sys_hubbard_t(hubbard)

        ! Deallocate all allocated components of a sys_hubbard_t object

        ! In/Out:
        !   sys_hubbard: sys_hubbard_t object to be deallocated.

        use system, only: sys_hubbard_t
        use symmetry_types, only: dealloc_mom_sym_t

        type(sys_hubbard_t), intent(inout) :: hubbard

        call dealloc_mom_sym_t(hubbard%mom_sym)

    end subroutine dealloc_sys_hubbard_t

    subroutine dealloc_sys_heisenberg_t(heisenberg)

        ! Deallocate all allocated components of a sys_heisenberg_t object

        ! In/Out:
        !   sys_heisenberg: sys_heisenberg_t object to be deallocated.

        use system, only: sys_heisenberg_t

        type(sys_heisenberg_t), intent(inout) :: heisenberg

        if (allocated(heisenberg%lattice_mask)) deallocate(heisenberg%lattice_mask)

    end subroutine dealloc_sys_heisenberg_t

    subroutine dealloc_sys_ueg_t(ueg)

        ! Deallocate all allocated components of a sys_ueg_t object

        ! In/Out:
        !   sys_ueg: sys_ueg_t object to be deallocated.

        use system, only: sys_ueg_t

       type(sys_ueg_t), intent(inout) :: ueg

        nullify(ueg%coulomb_int)
        nullify(ueg%exchange_int)

    end subroutine dealloc_sys_ueg_t

    subroutine dealloc_sys_read_in_t(read_in)

        ! Deallocate all allocated components of a sys_read_in_t object

        ! In/Out:
        !   sys_read_in: sys_read_in_t object to be deallocated.

        use system, only: sys_read_in_t
        use molecular_integrals, only: end_one_body_t, end_two_body_t
        use symmetry_types, only: dealloc_pg_sym_t

        type(sys_read_in_t), intent(inout) :: read_in

        call end_one_body_t(read_in%one_e_h_integrals)
        call end_one_body_t(read_in%one_body_op_integrals)
        call end_two_body_t(read_in%coulomb_integrals)
        call dealloc_pg_sym_t(read_in%pg_sym)

    end subroutine dealloc_sys_read_in_t

end module dealloc
