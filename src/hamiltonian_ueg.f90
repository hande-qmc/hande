module hamiltonian_ueg

! Module for evaluating Hamiltonian matrix elements for the uniform electron
! gas.

use const

implicit none

contains

    pure function get_hmatel_ueg(sys, f1, f2) result(hmatel)

        ! In:
        !    sys: system to be studied.
        !    f1, f2: bit string representation of the Slater
        !        determinants D1 and D2 respectively.
        ! Returns:
        !    Hamiltonian matrix element between the two determinants,
        !    < D1 | H | D2 >, where the determinants are formed from
        !    real space basis functions.

        ! Used in the UEG only.

        use excitations, only: excit_t, get_excitation
        use system, only: sys_t
        use hamiltonian_data

        type(hmatel_t) :: hmatel
        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f1(sys%basis%tot_string_len), f2(sys%basis%tot_string_len)
        type(excit_t) :: excitation

        ! Test to see if Hamiltonian matrix element is non-zero.

        ! Assume D1 and D2 are of the same symmetry.  Namely:

        !     We assume Ms is conserved (ie has already been checked for).

        !     The overall crystal momentum must be conserved (i.e. satisfy
        !     translational symmetry).  We assume this is also already checked.

        excitation = get_excitation(sys%nel, sys%basis, f1,f2)

        ! Connected determinants can differ by (at most) 2 spin orbitals.
        ! UEG (at least in the RHF basis of plane waves) has only double
        ! excitations, as the kinetic operator is diagonal in a plane wave
        ! basis.
        select case(excitation%nexcit)
        ! Apply Slater--Condon rules.
        case(0)

            ! < D | H | D > = \sum_i < i | h(i) | i > + \sum_i \sum_{j>i} < ij || ij >
            hmatel%r = slater_condon0_ueg(sys, f1)

        case(2)

            ! < D | H | D_{ij}^{ab} > = < ij || ab >

            ! Two electron operator
            hmatel%r = slater_condon2_ueg(sys, excitation%from_orb(1), excitation%from_orb(2), &
                                      & excitation%to_orb(1), excitation%to_orb(2),excitation%perm)
        case default

            hmatel%r = 0.0_p

        end select

    end function get_hmatel_ueg

    pure function slater_condon0_ueg(sys, f) result(hmatel)

        ! In:
        !    sys: system to be studied.
        !    f: bit string representation of the Slater determinant.
        ! Returns:
        !    < D_i | H | D_i >, the diagonal Hamiltonian matrix elements, for
        !        the Hubbard model in momentum space.

        use determinants, only: decode_det, sum_sp_eigenvalues_occ_list
        use system, only: sys_t

        real(p) :: hmatel
        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%tot_string_len)
        integer :: occ_list(sys%nel)

        call decode_det(sys%basis, f, occ_list)

        ! < D | H | D > = \sum_i < i | h(i) | i > + \sum_i \sum_{j>i} < ij || ij >
        hmatel = 0.0_p

        ! One electron operator: kinetic term
        hmatel = sum_sp_eigenvalues_occ_list(sys, occ_list)

        ! Two electron operator: Coulomb term.
        hmatel = hmatel + exchange_energy_ueg(sys, occ_list)

    end function slater_condon0_ueg

    pure function kinetic_energy_ueg(sys, f) result(hmatel)

        ! Calculate the kinetic energy of a given determinant.

        ! In:
        !    sys: system to be studied.
        !    f: bit string representation of the Slater determinant.
        ! Returns:
        !    < D_i | T | D_i >, the kinetic energy for the ueg.

        use determinants, only: decode_det, sum_sp_eigenvalues_occ_list
        use system, only: sys_t

        real(p) :: hmatel
        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%tot_string_len)
        integer :: occ_list(sys%nel)

        call decode_det(sys%basis, f, occ_list)

        ! < D | T | D > = \sum_i < i | h(i) | i >
        ! One electron operator: kinetic term
        hmatel = sum_sp_eigenvalues_occ_list(sys, occ_list)

    end function kinetic_energy_ueg

    pure function exchange_energy_ueg(sys, occ_list) result(hmatel)

        ! Cacluate the exchange energy from an orbital list.

        ! In:
        !    sys: system being studied.
        !    occ_list: list of occupied orbitals.

        use system, only: sys_t

        integer :: i, j
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: occ_list(:)
        real(p) :: hmatel

        hmatel = 0.0_p

        do i = 1, sys%nel
            do j = i+1, sys%nel
                ! Coulomb term is infinite but cancels exactly with the
                ! infinities in the electron-background and
                ! background-background interactions.
                if (mod(occ_list(i),2) == mod(occ_list(j),2)) then
                    ! Have an exchange term
                    hmatel = hmatel - sys%ueg%exchange_int(sys%lattice%box_length(1), sys%basis, occ_list(i), occ_list(j))
                end if
            end do
        end do

    end function exchange_energy_ueg

    pure function slater_condon2_ueg(sys, i, j, a, b, perm) result(hmatel)

        ! In:
        !    sys: system to be studied.
        !    i,j:  index of the spin-orbital from which an electron is excited in
        !          the reference determinant.
        !    a,b:  index of the spin-orbital into which an electron is excited in
        !          the excited determinant.
        !    perm: true if D and D_i^a are connected by an odd number of
        !          permutations.
        ! Returns:
        !    < D | H | D_ij^ab >, the Hamiltonian matrix element between a
        !    determinant and a double excitation of it in the UEG.

        use ueg_system, only: get_two_e_int_ueg
        use system, only: sys_t

        real(p) :: hmatel
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: i, j, a, b
        logical, intent(in) :: perm

        hmatel = get_two_e_int_ueg(sys, i, j, a, b)

        if (perm) hmatel = -hmatel

    end function slater_condon2_ueg

    pure function slater_condon2_ueg_excit(sys, i, a, b, perm) result(hmatel)

        ! In:
        !    sys: system to be studied.
        !    i:  index of the spin-orbital from which an electron is excited in
        !          the reference determinant.
        !    a,b:  index of the spin-orbital into which an electron is excited in
        !          the excited determinant.
        !    perm: true if D and D_i^a are connected by an odd number of
        !          permutations.
        ! Returns:
        !    < D | H | D_ij^ab >, the Hamiltonian matrix element between a
        !    determinant and a double excitation of it in the UEG, where
        !    j is defined such that momentum is conserved.

        ! WARNING: This function assumes that the D_{ij}^{ab} is a symmetry allowed
        ! excitation from D (and so the matrix element is *not* zero by
        ! symmetry).  This is less safe that slater_condon2_ueg but much faster
        ! as it allows symmetry checking to be skipped in the integral
        ! calculation.

        use system, only: sys_t

        real(p) :: hmatel
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: i, a, b
        logical, intent(in) :: perm

        hmatel = 0.0_p

        if (sys%basis%basis_fns(i)%Ms == sys%basis%basis_fns(a)%Ms) &
            hmatel = sys%ueg%coulomb_int(sys%lattice%box_length(1), sys%basis, i, a)
        if (sys%basis%basis_fns(i)%Ms == sys%basis%basis_fns(b)%Ms) &
            hmatel = hmatel - sys%ueg%coulomb_int(sys%lattice%box_length(1), sys%basis, i, b)

        if (perm) hmatel = -hmatel

    end function slater_condon2_ueg_excit

    subroutine create_weighted_excitation_list_ueg(sys, i, a_list, a_list_len, weights, weight_tot)

        ! Generate a list of allowed excitations from a to an element in a_list with their weights based on
        ! |<ia|ai>|. 
        !
        ! In:
        !    sys:   The system in which the orbitals live
        !    i:  integer specifying the from orbital
        !    a_list:   list of integers specifying the basis functions we're allowed to excite to
        !    a_list_len:   The length of a_list
        ! Out:
        !    weights:   A list of reals (length a_list_len), with the weight of each of the to_list orbitals
        !    weight_tot: The sum of all the weights.
        
        use system, only: sys_t
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: i, a_list_len, a_list(a_list_len)
        real(p), intent(out) :: weights(a_list_len), weight_tot
        
        integer :: k
        real(p) :: weight

        weight_tot = 0.0_p 
        do k=1, a_list_len
            if (a_list(k) /= i) then
                weight = sys%ueg%coulomb_int(sys%lattice%box_length(1), sys%basis, i, a_list(k))
                weight = abs(weight)
            else
                weight = 0.0_p
            end if
            weights(k) = weight
            weight_tot = weight_tot + weight
        end do

    end subroutine create_weighted_excitation_list_ueg

    pure function potential_energy_ueg(sys, f1, f2, excitation) result (potential_energy)

        ! In:
        !    sys: system of interest.
        !    f1, f2: bit string representations of two determinants.
        !    excitation: excit_t object describing the excitation connecting f1 and f2.
        ! Returns:
        !    The potential energy matrix element from a given bra (f1) and
        !    ket(f2), i.e. if H = T + V, then <f2| V |f1>.  This amounts to
        !    calculating the exchange contribution if f1=f2 and <D| H | D^ab>
        !    if f1 and f2 are related by a double excitation.

        use system, only: sys_t
        use excitations, only: excit_t
        use determinants, only: decode_det

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f1(:), f2(:)
        type(excit_t), intent(in) :: excitation

        real(p) :: potential_energy
        integer :: occ_list(sys%nel)

        potential_energy = 0.0_p

        select case(excitation%nexcit)
        case(0)
            ! Evaluate the exchange contribution.
            call decode_det(sys%basis, f1, occ_list)
            potential_energy = exchange_energy_ueg(sys, occ_list)
        case(2)
            ! Evaluate  <D| H | D^ab>.
            potential_energy = slater_condon2_ueg(sys, excitation%from_orb(1), excitation%from_orb(2), &
                                                  & excitation%to_orb(1), excitation%to_orb(2),excitation%perm)
        end select

    end function potential_energy_ueg

    function exchange_energy_orb(sys, occ_list, i) result (ex)

        ! Work out exchange energy for one orbital in particular.

        ! In:
        !    sys: system being studied.
        !    occ_list: list of occupied orbitals.
        !    i: orbital under consideration.
        ! Returns: Coulomb exchange integral ~ \sum_{iorb!=i} 1/|k_iorb - k_i|^2.


        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: occ_list(:)
        integer, intent(in) :: i

        real(p) :: ex
        integer :: iorb

        ex = 0.0_p
        do iorb = 1, sys%nel
            if (mod(occ_list(iorb),2) == mod(i,2) .and. occ_list(iorb) /= i) then
                ex = ex - sys%ueg%exchange_int(sys%lattice%box_length(1), sys%basis, occ_list(iorb), i)
            end if
        end do

    end function exchange_energy_orb

    function madelung_orb(sys, ref_occ_list, i) result (mad)

        ! Madelung contribution per electron.
        ! Expression fitted by Schoof et al., Phys. Rev. Lett. 115, 130402 (2015) using Fraser et al. Phys. Rev. B 53, 1814 (1996).

        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: ref_occ_list(:)
        integer, intent(in) :: i

        real(p) :: mad
        logical :: in_ref
        integer :: j

        ! Check whether orbital i is in reference.
        in_ref = .false.
        do j = 1, sys%nel
            if (ref_occ_list(j) == i) then
                in_ref = .true.
            end if
        end do

        if (in_ref) then
            mad = -2.837297*(0.75_p/(pi * (sys%ueg%r_s**3) * real(sys%nel,p)))**(1.0_p/3.0_p)
        else
            mad = 0.0_p
        end if


    end function madelung_orb

    subroutine calc_fock_values_3d_ueg(sys, propagator, ref_occ_list)
        
        use system, only: sys_t
        use qmc_data, only: propagator_t

        type(sys_t), intent(in) :: sys
        type(propagator_t), intent(inout) :: propagator
        integer, intent(in) :: ref_occ_list(:)
        
        integer :: iorb
        
        ! [todo]: Probably does not work properly if there is CAS since some ref det orbs are frozen!
        do iorb = 1, sys%basis%nbasis
            propagator%sp_fock(iorb) = propagator%sp_fock(iorb) + exchange_energy_orb(sys, ref_occ_list, iorb) + &
                0.5_p*madelung_orb(sys, ref_occ_list, iorb)
        end do

    end subroutine calc_fock_values_3d_ueg

end module hamiltonian_ueg
