module heisenberg_estimators

use const

implicit none

contains

    pure subroutine update_proj_energy_heisenberg_basic(sys, f0, cdet, pop, D0_pop_sum, proj_energy_sum, excitation, hmatel)

        ! Add the contribution of the current basis function to the
        ! projected energy.
        ! The projected energy estimator is:
        !   \sum_{i \neq 0} <D_i|H|D_0> N_i/N_0
        ! where N_i is the population on the i-th basis_function, D_i,
        ! and 0 refers to the reference basis function.
        ! During a MC cycle we store N_0 and \sum_{i \neq 0} <D_i|H|D_0> N_i
        ! This procedure is for the Heisenberg model only.

        ! In:
        !    sys: system being studied.
        !    f0: reference basis function.
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.  Only the bit string field needs to be set.
        !    pop: population on current determinant.
        ! In/Out:
        !    D0_pop_sum: running total of N_0, the population on the reference
        !        determinant, |D_0>.  Updated only if cdet is |D_0>.
        !    proj_energy_sum: running total of \sum_{i \neq 0} <D_i|H|D_0> N_i.
        !        Updated only if <D_i|H|D_0> is non-zero.
        !    excitation: excitation connecting the spin product to the trial wavefunction.
        ! Out:
        !    hmatel: <D_i|H|D_0>, the Hamiltonian matrix element between the
        !       spin product and the trial wavefunction.

        ! NOTE: it is the programmer's responsibility to ensure D0_pop_sum and
        ! proj_energy_sum are zero before the first call.

        use determinants, only: det_info_t
        use excitations, only: excit_t
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f0(:)
        type(det_info_t), intent(in) :: cdet
        real(p), intent(in) :: pop
        real(p), intent(inout) :: D0_pop_sum, proj_energy_sum
        type(excit_t), intent(inout) :: excitation
        real(p), intent(out) :: hmatel

        integer :: bit_position, bit_element

        hmatel = 0.0_p

        if (excitation%nexcit == 0) then
            ! Have reference determinant.
            D0_pop_sum = D0_pop_sum + pop
        else if (excitation%nexcit == 1) then
            ! Have a determinant connected to the reference determinant: add to
            ! projected energy.

            bit_position = sys%basis%bit_lookup(1,excitation%from_orb(1))
            bit_element = sys%basis%bit_lookup(2,excitation%from_orb(1))

            if (btest(sys%real_lattice%connected_orbs(bit_element, excitation%to_orb(1)), bit_position)) then
                 hmatel = -sys%heisenberg%J/2
                 proj_energy_sum = proj_energy_sum + hmatel*pop
             end if
        end if

    end subroutine update_proj_energy_heisenberg_basic

    pure subroutine update_proj_energy_heisenberg_neel_singlet(sys, f0, cdet, pop, D0_pop_sum, proj_energy_sum, excitation, hmatel)

        ! Add the contribution of the current basis function to the
        ! projected energy.
        ! This uses the Neel singlet state as a trial function.
        ! This function is an integral over the Neel states in all directions
        ! (giving equal weight to each). Hence it is rotationally invariant
        ! and hence is an S=0 eigenstate. The ground state is known to be an
        ! S = 0 eigenstate also (ie of the total spin squared operator). This
        ! makes this state a suitable trial function. For details on the state, see
        ! K. Runge, Phys. Rev. B 45, 7229 (1992).
        ! This procedure is for the Heisenberg model only.

        ! In:
        !    sys: system being studied.
        !    f0: reference basis function (unused, for interface compatibility only).
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.  Only the bit string and data fields need to be set.
        !    pop: population on current determinant, unweighted by the trial function.
        ! In/Out:
        !    D0_pop_sum: running total of N_0, the population on the reference
        !        determinant, |D_0>.  Updated only if cdet is |D_0>.
        !    proj_energy_sum: running total of \sum_{i \neq 0} <D_i|H|D_0> N_i.
        !        Updated only if <D_i|H|D_0> is non-zero.
        !    excitation: excitation connecting the spin product to the trial wavefunction.
        !       As each spin product is in the trial wavefunction, this is
        !       simply null, but included for interface compatibility.
        ! Out:
        !    hmatel: <D_i|H|D_0>, the Hamiltonian matrix element between the
        !       spin product and the trial wavefunction.

        ! NOTE:

        ! 1. It is the programmer's responsibility to ensure D0_pop_sum and
        !    proj_energy_sum are zero before the first call.

        ! 2. With importance sampling, the psip amplitudes, n_i, represent the quantities
        !    f_i = c_i^T * c_i
        !    where c_i the amplitude of |D_i> in the true ground state and c_i^T
        !    is the amplitude of |D_i> in the trial ground state.   pop represents
        !    c_i rather than f_i.

        use determinants, only: det_info_t
        use excitations, only: excit_t
        use fciqmc_data, only: neel_singlet_amp
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f0(:)
        type(det_info_t), intent(in) :: cdet
        real(p), intent(in) :: pop
        real(p), intent(inout) :: D0_pop_sum, proj_energy_sum
        type(excit_t), intent(inout) :: excitation
        real(p), intent(out) :: hmatel

        integer :: n, lattice_1_up, lattice_2_up

        excitation = excit_t(0, (/ 0,0 /), (/ 0,0 /), .false.)

        ! WARNING: assuming trial function parameters are the last two elements
        ! in cdet%data.
        n = nint(cdet%data(size(cdet%data)-1))
        lattice_1_up = nint(cdet%data(size(cdet%data)))

        ! Deduce the number of 0-1 bonds where the 1's are on the second sublattice:
        ! The total number of 0-1 bonds, n(0-1), can be found from the diagonal
        ! element of the current basis function:
        ! hmatel = -J*(nbonds - 2*n(0-1))/4
        ! This means we can avoid calculating n(0-1) again, which is expensive.
        ! We know the number of 0-1 bonds where the 1 (the spin up) is on sublattice 1,
        ! so can then deduce the number where the 1 is on sublattice 2.
        lattice_2_up = ((sys%heisenberg%nbonds) + nint(4*cdet%data(1)/sys%heisenberg%J))/2 - lattice_1_up

        ! There are three contributions to add to the projected energy from
        ! the current basis function. Consider the Neel singlet state:
        ! |NS> = \sum_{i} a_i|D_i>
        ! The amplitudes, a_i, only depend on the number of spins up on sublattice 1.
        ! We want to calculate \sum_{i} (a_i * <D_i|H|D_j> * n_j), where |D_j> is the
        ! current basis function, and then add this to the current projected energy.

        ! Firstly, consider the diagonal term:
        ! We have <D_j|H|D_j> stored, so this is simple:
        hmatel = neel_singlet_amp(n) * cdet%data(1)

        ! Now consider off-diagonal contributions, |D_i> /= |D_j>.
        ! To create a connected |D_i> from |D_j>, simply find a pair of anti-parallel
        ! spins in |D_j> (a 0-1 bond) and flip both of these spins. Following this
        ! procedure for every pair of anti-parallel spins generates *all* connected
        ! |D_i>'s, and hence every term which needs to be considered for the following.
        ! Now, the amplitude a_i only depends on the number of spins up on sublattice 1.
        ! This will depend on whether, for the 0-1 bond flipped, the 1 (up spin) was one
        ! sublattice 1 or 2. If the up spin was on lattice 1, there will be one less up
        ! spin on sublattice 1 after the flip. If it was on sublattice 2, there will be
        ! one extra spin. So we can have either of the two amplitudes, neel_singlet_amp(n-1)
        ! or neel_singlet_amp(n+1) respectively. We just need to know how many of each type
        ! of 0-1 bonds there are. But we have these already - they are stored as lattice_1_up
        ! and lattice_2_up. Finally note that the matrix element is -J/2, and we can put
        ! this together...

        ! From 0-1 bonds where the 1 is on sublattice 1, we have:
        hmatel = hmatel - (sys%heisenberg%J * lattice_1_up * neel_singlet_amp(n-1))/2

        ! And from 1-0 bond where the 1 is on sublattice 2, we have:
        hmatel = hmatel - (sys%heisenberg%J * lattice_2_up * neel_singlet_amp(n+1))/2

        proj_energy_sum = proj_energy_sum + hmatel * pop

        ! Now we just need to find the contribution to the denominator. The total
        ! denominator is
        ! \sum_{i} (a_i * n_i)
        ! Hence from this particular basis function, |D_j>, we just add (a_j * n_j)

        D0_pop_sum = D0_pop_sum + pop*neel_singlet_amp(n)

    end subroutine update_proj_energy_heisenberg_neel_singlet

    pure function neel_singlet_data(sys, f) result(spin_config_data)

        ! This subroutine calculates the number of up spins on one of the
        ! sublattices so that it can be stored, and also the total number
        ! of 0-1 bonds where the up spins lie on this same sublattice. These are
        ! later used in the update_proj_energy_heisenberg_neel_singlet subroutine.
        ! This procedure is for the Heisenberg model only.

        ! In:
        !    sys: system being studied.
        !    f: bit string representation of the basis function.
        ! Returns:
        !    spin_config_data: A two component integer vector, where the first
        !       component stores the total number of spins up on the first sublattice
        !       for this basis function, and the second component stores the numbers
        !       of 0-1 bonds, where the 1 (the up spin) is on the first sublattice.

        use bit_utils, only: count_set_bits
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%string_len)
        integer(i0) :: f_mask(sys%basis%string_len), f_not(sys%basis%string_len), g(sys%basis%string_len)
        integer :: lattice_1_up, n, i, ipos, basis_find
        integer :: spin_config_data(2)

        ! Calculate the number of up spins on the first sublattice.
        f_mask = iand(f, sys%heisenberg%lattice_mask)
        n = sum(count_set_bits(f_mask))

        ! Find the number of 0-1 bonds where the 1 lies on the first sublattice.
        lattice_1_up = 0
        f_not = not(f)
        do i = 1, sys%basis%string_len
            do ipos = 0, i0_end
                if (btest(f_mask(i), ipos)) then
                    basis_find = sys%basis%basis_lookup(ipos, i)
                    g = iand(f_not, sys%real_lattice%connected_orbs(:,basis_find))
                    lattice_1_up = lattice_1_up + sum(count_set_bits(g))
                end if
            end do
        end do

        spin_config_data(1) = n
        spin_config_data(2) = lattice_1_up

    end function neel_singlet_data

end module heisenberg_estimators
