module hamiltonian

use const
use basis
use determinants
use parallel

implicit none

contains

    elemental function get_hmatel(d1, d2) result(hmatel)

        ! In:
        !    d1, d2: integer labels of two determinants, as stored in the
        !            dets array.
        ! Returns:
        !    Hamiltonian matrix element between the two determinants, 
        !    < D1 | H | D2 >.

        ! This is just a wrapper function around the system specific get_hmatel
        ! functions.

        ! Having separate functions for the different systems might seem
        ! somewhat redundant (a lot of the code in the functions is similar)
        ! but enables us to use only one test for the system type.  A small
        ! efficiency for not much effort. :-)

        real(dp) :: hmatel
        integer, intent(in) :: d1, d2

        select case(system_type)
        case(hub_k)
            hmatel = get_hmatel_k(d1, d2)
        case(hub_real)
            hmatel = get_hmatel_real(d1, d2)
        end select

    end function get_hmatel

    elemental function get_hmatel_k(d1, d2) result(hmatel)

        ! In:
        !    d1, d2: integer labels of two determinants, as stored in the
        !            dets_* arrays.
        ! Returns:
        !    Hamiltonian matrix element between the two determinants, 
        !    < D1 | H | D2 >, where the determinants are formed from
        !    momentum space basis functions.

        ! Used in the momentum space formulation of the Hubbard model only.

        use hubbard_k

        real(dp) :: hmatel
        integer, intent(in) :: d1, d2
        logical :: non_zero
        type(excit) :: excitation
        integer :: root_det(nel)
        integer :: i, j

        hmatel = 0.0_dp
        non_zero = .false.

        ! Test to see if Hamiltonian matrix element is non-zero.

        ! Assume D1 and D2 are of the same symmetry.  Namely:

            ! We assume Ms is conserved (ie has already been checked for).

            ! In the momentum space description the overall crystal 
            ! momentum must be conserved up to a reciprocal lattice
            ! vector (i.e. satisfy translational symmetry).
            ! We assume this is also already checked.

        excitation = get_excitation(dets_list(:,d1), dets_list(:,d2))
        ! Connected determinants can differ by (at most) 2 spin orbitals.
        if (excitation%nexcit <= 2) then
            non_zero = .true.
        end if

        if (non_zero) then
            select case(excitation%nexcit)
            ! Apply Slater--Condon rules.
            case(0)

                root_det = decode_det(dets_list(:,d1))

                ! < D | H | D > = \sum_i < i | h(i) | i > + \sum_i \sum_{j>i} < ij || ij >

                ! One electron operator
                do i = 1, nel
                    hmatel = hmatel + get_one_e_int_k(root_det(i), root_det(i))
                end do

                ! Two electron operator
                do i = 1, nel
                    do j = i+1, nel
                        hmatel = hmatel + get_two_e_int_k(root_det(i), root_det(j), root_det(i), root_det(j))
                    end do
                end do

!            case(1)

                ! < D | H | D_i^a > = < i | h(a) | a > + \sum_j < ij || aj >

                ! Single excitations are not connected in the momentum space
                ! basis.

                ! One electron operator
                ! The kinetic operator is diagonal in the momentum space basis.

                ! Two electron operator
                ! < ij | aj > = 0 only if crystal momentum is conserved up to
                ! a reciprocal lattice vector.
                ! As k_i /= k_j, this cannot be met.

            case(2)

                ! < D | H | D_{ij}^{ab} > = < ij || ab >

                ! Two electron operator
                hmatel = get_two_e_int_k(excitation%from_orb(1), excitation%from_orb(2), excitation%to_orb(1), excitation%to_orb(2))

                if (excitation%perm) hmatel = -hmatel

            end select
        end if

    end function get_hmatel_k

    elemental function get_hmatel_real(d1, d2) result(hmatel)

        ! In:
        !    d1, d2: integer labels of two determinants, as stored in the
        !            dets array.
        ! Returns:
        !    Hamiltonian matrix element between the two determinants, 
        !    < D1 | H | D2 >, where the determinants are formed from
        !    real space basis functions.

        ! Used in the real space formulation of the Hubbard model only.

        use hubbard_real

        real(dp) :: hmatel
        integer, intent(in) :: d1, d2
        logical :: non_zero
        type(excit) :: excitation
        integer :: root_det(nel)
        integer :: i

        hmatel = 0.0_dp
        non_zero = .false.

        ! Test to see if Hamiltonian matrix element is non-zero.

        ! We assume Ms is conserved (ie has already been checked for).
        excitation = get_excitation(dets_list(:,d1), dets_list(:,d2))
        ! Connected determinants can differ by (at most) 2 spin orbitals.
        if (excitation%nexcit <= 2) then
            ! Space group symmetry not currently implemented.
            non_zero = .true.
        end if

        ! Matrix elements in the real space formulation are quite simple.

        ! 1. < i | T | i > = 0
        !    Thus the one-electron terms only occur between single excitation
        !    matrix elements.
        ! 2. < m,s1 n,s2 | U | p,s1 q,s2 > = U \delta_{m,n} \delta_{m,p} \delta_{m,q} , s1/=s2
        !    Thus the Coulomb integrals that occur in < D | H | D_i^a > and
        !    < D | H | D_{ij}^{ab} > are zero.

        if (non_zero) then
            select case(excitation%nexcit)
            ! Apply Slater--Condon rules.
            case(0)

                ! < D | H | D > = \sum_i < i | h(i) | i > + \sum_i \sum_{j>i} < ij || ij >

                ! < i | T | i > = 0 within the real space formulation of the
                ! Hubbard model, unless site i is its own periodic image, in
                ! which case it has a kinetic interaction with its self-image.
                ! This only arises if there is at least one crystal cell vector
                ! which is a unit cell vector.
                if (t_self_images) then
                    root_det = decode_det(dets_list(:,d1))
                    do i = 1, nel
                        hmatel = hmatel + get_one_e_int_real(root_det(i), root_det(i))
                    end do
                end if

                ! Two electron operator
                hmatel = hmatel + get_coulomb_matel_real(dets_list(:,d1))

            case(1)

                ! < D | H | D_i^a > = < i | h(a) | a > + \sum_j < ij || aj >

                ! One electron operator
                 hmatel = get_one_e_int_real(excitation%from_orb(1), excitation%to_orb(1)) 

                ! Two electron operator
                ! < D | U | D_i^a > = 0 within the real space formulation of the
                ! Hubbard model.

                if (excitation%perm) hmatel = -hmatel

!            case(2)

                ! < D | H | D_{ij}^{ab} > = < ij || ab >
                !                         = 0 within the real space formulation
                !                             of the Hubbard model.

            end select
        end if

    end function get_hmatel_real

end module hamiltonian
