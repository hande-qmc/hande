module canonical_kinetic_energy

! Estimate the canonical free-electron total energy.

use const

implicit none

! Number of report loops use to estimate the energy.
integer :: nkinetic_cycles

contains

    subroutine estimate_kinetic_energy(sys)

        ! From the Fermi factors calculated in the grand canonical ensemble we can
        ! estimate the total energy in the canonical ensemble by generating determinants
        ! with arbitrary particle number and only keeping those with <N> = nel.

        ! In:
        !    sys: system being studied.

        use system
        use dSFMT_interface, only: dSFMT_t, dSFMT_init
        use parallel
        use calc, only: ms_in, seed
        use fciqmc_data, only: init_beta, D0_population
        use utils, only: rng_init_info
        use hamiltonian_ueg, only: sum_sp_eigenvalues

        type(sys_t), intent(inout) :: sys

        real(dp) :: p_single(sys%basis%nbasis/2)
        real(dp) :: r
        integer :: occ_list(sys%nel)
        logical :: gen
        real(p) :: ke_proc, ke
        integer :: ierr, ireport, iaccept, iorb

        type(sys_t) :: sys_bak
        type (dSFMT_t) :: rng

        if (parent) call rng_init_info(seed+iproc)
        call dSFMT_init(seed+iproc, 50000, rng)
        call copy_sys_spin_info(sys, sys_bak)
        call set_spin_polarisation(sys%basis%nbasis, ms_in, sys)

        if (parent) write (6,'(1X,a12,3X,a19)') '# iterations', 'E_0'

        forall(iorb=1:sys%basis%nbasis:2) p_single(iorb/2+1) = 1.0_p / &
                                                          (1+exp(init_beta*(sys%basis%basis_fns(iorb)%sp_eigv-sys%ueg%chem_pot)))

        ! iaccept is analogous to mc_cycles in the normal qmc algorithms so we
        ! do nkinetic_cycles*init_pop iterations in total.
        do ireport = 1, nkinetic_cycles
            ke = 0.0_p
            ke_proc = 0.0_p
            iaccept = 0
            do while (iaccept < D0_population)
                if (sys%nalpha > 0) call generate_allowed_orbital_list(sys, rng, p_single, sys%nalpha, &
                                                                       1, occ_list(:sys%nalpha), gen)
                if (.not. gen) cycle
                if (sys%nbeta > 0) call generate_allowed_orbital_list(sys, rng, p_single, sys%nbeta, &
                                                                      0, occ_list(sys%nalpha+1:), gen)
                if (.not. gen) cycle
                iaccept = iaccept + 1
                ke_proc = ke_proc + sum_sp_eigenvalues(sys, occ_list)
            end do

#ifdef PARALLEL
            call mpi_reduce(ke_proc, ke, 1, mpi_preal, MPI_SUM, root, MPI_COMM_WORLD, ierr)
#else
            ke = ke_proc
#endif
            if (parent) write(6,'(3X,i10,5X,es17.10)') ireport, ke / (nprocs*D0_population)
        end do

        if (parent) write(6, '()')

    end subroutine estimate_kinetic_energy

    subroutine generate_allowed_orbital_list(sys, rng, porb, nselect, spin_factor, occ_list, gen)

        ! Generate a list of orbitals according to their single
        ! particle GC orbital occupancy probabilities.

        ! In:
        !    porb: porb(i) gives the probabilty of selecting
        !        the orbital i.
        !    nselect: number of orbitals to select.
        !    spin_factor: integer to account for odd/even ordering of
        !        alpha/beta spin orbitals. Set to 1 for alpha spins, 0 for beta spins.
        ! In/Out:
        !    rng: random number generator.
        !    occ_list: array containing occupied orbitals.

        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        real(dp), intent(in) :: porb(:)
        integer, intent(in) :: nselect
        integer, intent(in) :: spin_factor
        type(dSFMT_t), intent(inout) :: rng
        integer, intent(out) :: occ_list(:)
        logical, intent(out) :: gen

        integer :: iorb, iselect
        real(dp) :: r

        iselect = 0
        occ_list = 0

        do iorb = 1, sys%basis%nbasis/2
            ! Select a random orbital.
            r = get_rand_close_open(rng)
            if (porb(iorb) > r) then
                iselect = iselect + 1
                if (iselect > nselect) then
                    ! Selected too many.
                    gen = .false.
                    exit
                end if
                occ_list(iselect) = 2*iorb - spin_factor
            end if
        end do
        if (iselect == nselect) then
            gen = .true.
        else
            gen = .false.
        end if

    end subroutine generate_allowed_orbital_list

end module canonical_kinetic_energy
