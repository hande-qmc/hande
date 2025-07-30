module excit_gen_trotter

! Module for random excitation generators and related routines for the VTV-type
! Trotter error norm calculations, for systems read from an FCIDUMP file.

! These functions are based on those for general molecular Hamiltonians in
! excit_gen_mol.F90. However, the VTV commutator only contains terms that
! cause single excitations, and are three-body in nature. We define new
! functions here to allow possible optimizations to be made as needed.

use const, only: i0, p

implicit none

contains

    subroutine gen_excit_trotter(rng, sys, excit_gen_data, cdet, pgen, connection, hmatel, allowed_excitation)

        ! Create a random excitation from cdet and calculate both the probability
        ! of selecting that excitation and the Hamiltonian matrix element.

        ! In:
        !    sys: system object being studied.
        !    excit_gen_data: Data for the excitation generator.
        !    parent_sign: sign of the population on the parent determinant (i.e.
        !        either a positive or negative integer).
        ! In/Out:
        !    rng: random number generator.
        !    cdet: info on the current determinant (cdet) that we will gen from.
        ! Out:
        !    pgen: probability of generating the excited determinant from cdet.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are gened.
        !    hmatel: < D | H | D' >, the Hamiltonian matrix element between a
        !        determinant and a connected determinant.
        !    allowed_excitation: false if a valid symmetry allowed excitation
        !        was not generated.

        use determinant_data, only: det_info_t
        use excitations, only: excit_t
        use excit_gens, only: excit_gen_data_t
        use system, only: sys_t
        use hamiltonian_data

        use dSFMT_interface, only: dSFMT_t, get_rand_close_open

        type(sys_t), intent(in) :: sys
        type(excit_gen_data_t), intent(in) :: excit_gen_data
        type(det_info_t), intent(inout) :: cdet
        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(out) :: pgen
        type(hmatel_t), intent(out) :: hmatel
        type(excit_t), intent(out) :: connection
        logical, intent(out) :: allowed_excitation

        ! We only have single excitations for this commutator.
        call gen_single_excit_trotter(rng, sys, cdet, pgen, connection, hmatel, &
                                allowed_excitation)

    end subroutine gen_excit_trotter

    subroutine gen_single_excit_trotter(rng, sys, cdet, pgen, connection, hmatel, allowed_excitation)

        ! Create a random single excitation from cdet and calculate both the
        ! probability of selecting that excitation and the Hamiltonian matrix
        ! element.

        ! In:
        !    sys: system object being studied.
        !    cdet: info on the current determinant (cdet) that we will gen
        !        from.
        ! In/Out:
        !    rng: random number generator.
        ! Out:
        !    pgen: probability of generating the excited determinant from cdet.
        !    connection: excitation connection between the current determinant
        !        and the child determinant.
        !    hmatel: < D | H | D' >, the Hamiltonian matrix element between a
        !       determinant and a connected determinant.
        !    allowed_excitation: false if an allowed excitation was not generated.

        use determinant_data, only: det_info_t
        use excitations, only: excit_t
        use excitations, only: find_excitation_permutation1
        use proc_pointers, only: slater_condon1_excit_ptr
        use system, only: sys_t
        use hamiltonian_data, only: hmatel_t
        use dSFMT_interface, only: dSFMT_t

        type(sys_t), intent(in) :: sys
        type(det_info_t), intent(in) :: cdet
        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(out) :: pgen
        type(hmatel_t), intent(out) :: hmatel
        type(excit_t), intent(out) :: connection
        logical, intent(out) :: allowed_excitation

        ! Select orbital to excite from and orbital to excite into.
        call choose_ia_trotter(rng, sys, sys%read_in%pg_sym%gamma_sym, cdet%f, cdet%occ_list, cdet%symunocc, &
                            connection%from_orb(1), connection%to_orb(1), allowed_excitation)
        ! This is a single excitation.
        connection%nexcit = 1

        if (allowed_excitation) then
            ! Calculate the probability of generating this excitation.
            pgen = calc_pgen_single_trotter(sys, sys%read_in%pg_sym%gamma_sym, cdet%occ_list, &
                                                            cdet%symunocc, connection%to_orb(1))

            ! Find the connecting matrix element.
            hmatel = slater_condon1_excit_ptr(sys, cdet%occ_list, connection%from_orb(1), &
                                          connection%to_orb(1), connection%perm)
        else
            ! This determinant has no single excitations at all.
            ! We simply return a null excitation.
            hmatel%c = cmplx(0.0_p, 0.0_p, p)
            hmatel%r = 0.0_p
            pgen = 1.0_p
        end if

    end subroutine gen_single_excit_trotter

    subroutine choose_ia_trotter(rng, sys, op_sym, f, occ_list, symunocc, i, a, allowed_excitation)

        ! Randomly choose a single excitation, i->a, of a determinant.

        ! In:
        !    sys: system object being studied.
        !    op_sym: symmetry of connecting operator.
        !    f: bit string representation of the Slater determinant from which
        !        an electron is excited.
        !    occ_list: integer list of occupied spin orbitals in the determinant.
        !        (min length: sys%nel)
        !    symunocc: number of unoccupied orbitals of each spin and
        !        irreducible representation. The same indexing scheme as
        !        nbasis_sym_spin (in the point_group_symmetry module) is used.
        ! In/Out:
        !    rng: random number generator.
        ! Out:
        !    i: orbital in determinant from which an electron is excited.
        !    a: previously unoccupied orbital into which an electron is excited.
        !        Not set if allowed_excitation is false.
        !    allowed_excitation: false if there are no possible single
        !        excitations from the determinant which conserve spin and spatial
        !        symmetry.

        use system, only: sys_t

        use dSFMT_interface, only: dSFMT_t, get_rand_close_open

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: op_sym
        integer(i0), intent(in) :: f(sys%basis%tot_string_len)
        integer, intent(in) :: occ_list(:), symunocc(:,sys%sym0_tot:)
        type(dSFMT_t), intent(inout) :: rng
        integer, intent(out) :: i, a
        logical, intent(out) :: allowed_excitation

        integer :: imsa, isyma, ind

        ! Does this determinant have any possible single excitations?
        allowed_excitation = .false.
        do i = 1, sys%nel
            imsa = (sys%basis%basis_fns(occ_list(i))%Ms+3)/2
            ! Get the required symmetry of the orbital bing excited to.
            isyma = sys%read_in%cross_product_sym_ptr(sys%read_in, &
                    sys%basis%basis_fns(occ_list(i))%sym, op_sym)
            if (symunocc(imsa, isyma) /= 0) then
                allowed_excitation = .true.
                exit
            end if
        end do

        if (allowed_excitation) then
            do
                ! Select an occupied orbital at random.
                i = occ_list(int(get_rand_close_open(rng)*sys%nel)+1)
                ! Conserve symmetry (spatial and spin) in selecting a.
                imsa = (sys%basis%basis_fns(i)%Ms+3)/2
                ! Get the required symmetry of the orbital bing excited to.
                isyma = sys%read_in%cross_product_sym_ptr(sys%read_in, sys%basis%basis_fns(i)%sym, op_sym)
                if (symunocc(imsa, isyma) /= 0) then
                        ! Found i.  Now find a...
                        ! It's cheaper to draw additional random numbers than
                        ! decode the full list of unoccupied orbitals,
                        ! especially as the number of basis functions is usually
                        ! much larger than the number of electrons.
                        do
                            ind = int(sys%read_in%pg_sym%nbasis_sym_spin(imsa,isyma)*get_rand_close_open(rng))+1
                            a = sys%read_in%pg_sym%sym_spin_basis_fns(ind,imsa,isyma)
                            if (.not.btest(f(sys%basis%bit_lookup(2,a)), sys%basis%bit_lookup(1,a))) exit
                        end do
                    exit
                end if
            end do

        end if

    end subroutine choose_ia_trotter

    pure function calc_pgen_single_trotter(sys, op_sym, occ_list, symunocc, a) result(pgen)

        ! In:
        !    sys: system object being studied.
        !    op_sym: symmetry of connecting operator.
        !    occ_list: integer list of occupied spin-orbitals in the determinant.
        !        (min length: sys%nel.)
        !    symunocc: number of unoccupied orbitals of each spin and
        !        irreducible representation.  The same indexing scheme as
        !        nbasis_sym_spin (in the point_group_symmetry module) is used.
        !    a: previously unoccupied orbital into which an electron is excited.
        ! Returns:
        !    The probability of generating the excitation i->a, where we select
        !    i uniformly from the list of occupied orbitals with allowed
        !    excitations (i.e. at least one single excitation exists involving
        !    those orbitals) and a is selected from the list of unoccupied
        !    orbitals which conserve spin and spatial symmetry *and* given that
        !    a single excitation has been selected.

        ! WARNING: We assume that the excitation is actually valid.
        ! This routine does *not* calculate the correct probability that
        ! a forbidden excitation (e.g. no possible single excitations exist) is
        ! generated.  The correct way to handle those excitations is to
        ! explicitly reject them. The generation probabilites of allowed
        ! excitations correctly take into account such rejected events.

        use system, only: sys_t

        real(p) :: pgen
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: op_sym
        integer, intent(in) :: occ_list(:), symunocc(:,sys%sym0_tot:), a

        integer :: imsi, isymi, imsa, isyma
        integer :: i, ni

        ! The generation probability, pgen, is:
        !   pgen = p_single p(i) p(a|i)
        !        = p_single 1/n_i 1/symunocc(ms_i, sym_i)
        ! where n_i is the number of electrons which have at least one possbile
        ! excitation. Note that the probability of selecting a single is 1 for
        ! this Hamiltonian, hence p_single = 1.

        ni = sys%nel
        do i = 1, sys%nel
            imsa = (sys%basis%basis_fns(occ_list(i))%Ms+3)/2
            ! Get the required symmetry of the orbital bing excited to.
            isyma = sys%read_in%cross_product_sym_ptr(sys%read_in, sys%basis%basis_fns(occ_list(i))%sym, op_sym)
            if (symunocc(imsa,isyma) == 0) ni = ni - 1
        end do

        imsi = (sys%basis%basis_fns(a)%Ms+3)/2
        isymi = sys%basis%basis_fns(a)%sym
        pgen = 1.0_p/(ni*symunocc(imsi,isymi))

    end function calc_pgen_single_trotter

end module excit_gen_trotter
