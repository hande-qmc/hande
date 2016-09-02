module ccmc_data

use const, only: i0, p

implicit none

! Work around Fortran's lack of arrays of pointers...
type bit_string_ptr
    integer(i0), pointer :: f(:)
end type bit_string_ptr

! Information about a cluster of excitors in addition to that stored in
! a det_info_t variable.
type cluster_t
    ! Pointers to the determinants formed by applying individual excitors to the
    ! reference determinant.
    type(bit_string_ptr), allocatable :: excitors(:) ! max: ex_level+2
    ! Number of excitors in cluster
    integer :: nexcitors
    ! Excitation level relative to the reference determinant of the determinant
    ! formed by applying the cluster of excitors to the reference determinant.
    integer :: excitation_level
    ! Overall amplitude of the cluster.
    real(p) :: amplitude
    ! < D_i | a_i D_0 >, where D_i is the determinant formed by applying the
    ! cluster of excitors to the reference determinant.  Equal to +1 or -1.
    integer :: cluster_to_det_sign
    ! probability of selecting this cluster at random
    real(p) :: pselect
end type cluster_t

type multispawn_stats_t
    integer :: nevents = 0
    integer :: nspawnings = 0
    integer :: nspawnings_max = 0
end type multispawn_stats_t

contains

    subroutine ms_stats_update(nspawnings, ms_stats)

        ! Update a multispawn_stats_t object.

        ! In:
        !    nspawnings: number of spawning attempts to be made from a single cluster selection.
        ! In/Out:
        !    ms_stats: On output, nevents is incremented and the nspawnings counts are increased accordingly.

        integer, intent(in) :: nspawnings
        type(multispawn_stats_t), intent(inout) :: ms_stats

        if (nspawnings > 1) then
            ms_stats%nevents = ms_stats%nevents + 1
            ms_stats%nspawnings = ms_stats%nspawnings + nspawnings
            ms_stats%nspawnings_max = max(ms_stats%nspawnings_max,nspawnings)
        end if

    end subroutine ms_stats_update

    pure function ms_stats_reduction(ms_stats_arr) result(ms_stats)

        ! In:
        !    ms_stats_arr: array of multispawn_stats_t objects.
        ! Returns:
        !    The multispawn_stats_t overall object (ie totals and max values from each component).

        type(multispawn_stats_t) :: ms_stats
        type(multispawn_stats_t), intent(in) :: ms_stats_arr(:)

        ms_stats%nevents = sum(ms_stats_arr%nevents)
        ms_stats%nspawnings = sum(ms_stats_arr%nspawnings)
        ms_stats%nspawnings_max = maxval(ms_stats_arr%nspawnings_max)

    end function ms_stats_reduction

    subroutine multispawn_stats_report(ms_stats)

        ! Print out a report of the cluster multispawn events.

        ! In:
        !    ms_stats: array of multispawn_stats_t objects.  The reduced version is printed out.

        use parallel
        use utils, only: int_fmt

        type(multispawn_stats_t), intent(in) :: ms_stats(:)
        type(multispawn_stats_t) :: ms_stats_total
#ifdef PARALLEL
        type(multispawn_stats_t) :: ms_stats_local
        integer :: ierr

        ms_stats_local = ms_stats_reduction(ms_stats)
        call mpi_reduce(ms_stats_local%nevents, ms_stats_total%nevents, 1, MPI_INTEGER, &
                        MPI_SUM, root, MPI_COMM_WORLD, ierr)
        call mpi_reduce(ms_stats_local%nspawnings, ms_stats_total%nspawnings, 1, MPI_INTEGER, &
                        MPI_SUM, root, MPI_COMM_WORLD, ierr)
        call mpi_reduce(ms_stats_local%nspawnings_max, ms_stats_total%nspawnings_max, 1, MPI_INTEGER, &
                        MPI_MAX, root, MPI_COMM_WORLD, ierr)
#else
        ms_stats_total = ms_stats_reduction(ms_stats)
#endif
        if (ms_stats_total%nevents > 0 .and. parent) then
            write (6,'(1X,"Multiple spawning events occurred.")')
            write (6,'(1X,"Number of multiple spawning events:",'//int_fmt(ms_stats_total%nevents,1)//')') &
                ms_stats_total%nevents
            write (6,'(1X,"Mean number of multiple spawning attempts per event:",2X,f11.2)') &
                real(ms_stats_total%nspawnings)/ms_stats_total%nevents
            write (6,'(1X,"Largest multiple spawning in a single event:",'//int_fmt(ms_stats_total%nspawnings_max,1)//',/)') &
                ms_stats_total%nspawnings_max
        end if

    end subroutine multispawn_stats_report

    pure subroutine convert_excitor_to_determinant(excitor, excitor_level, excitor_sign, f)

        ! We usually consider an excitor as a bit string representation of the
        ! determinant formed by applying the excitor (a group of annihilation
        ! and creation operators) to the reference determinant; indeed the
        ! determinant form is required when constructing Hamiltonian matrix
        ! elements.  However, the resulting determinant might actually contain
        ! a sign change, which needs to be absorbed into the (signed) population
        ! of excips on the excitor.

        ! This results from the fact that a determinant, |D>, is defined as:
        !   |D> = a^+_i a^+_j ... a^+_k |0>,
        ! where |0> is the vacuum, a^+_i creates an electron in the i-th
        ! orbital, i<j<...<k and |0> is the vacuum.  An excitor is defined as
        !   t_{ij...k}^{ab...c} = a^+_a a^+_b ... a^+_c a_k ... a_j a_i
        ! where i<j<...<k and a<b<...<c.  (This definition is somewhat
        ! arbitrary; the key thing is to be consistent.)  Hence applying an
        ! excitor to the reference might result in a change of sign, i.e.
        ! t_{ij}^{ab} |D_0> = -|D_{ij}^{ab}>.  As a more concrete example,
        ! consider a set of spin states (as the extension to fermions is
        ! irrelevant to the argument), with the reference:
        !   |D_0> = | 1 2 3 >
        ! and the excitor
        !   t_{13}^{58}
        ! Thus, using |0> to denote the vacuum:
        !   t_{13}^{58} | 1 2 3 > = + a^+_5 a^+_8 a_3 a_1 a^+_1 a^+_2 a^+_3 |0>
        !                         = + a^+_5 a^+_8 a_3 a^+_2 a^+_3 |0>
        !                         = - a^+_5 a^+_8 a_3 a^+_3 a^+_2 |0>
        !                         = - a^+_5 a^+_8 a^+_2 |0>
        !                         = + a^+_5 a^+_2 a^+_8 |0>
        !                         = - a^+_2 a^+_5 a^+_8 |0>
        !                         = - | 2 5 8 >
        ! Similarly
        !   t_{12}^{58} | 1 2 3 > = + a^+_5 a^+_8 a_2 a_1 a^+_1 a^+_2 a^+_3 |0>
        !                         = + a^+_5 a^+_8 a_2 a^+_2 a^+_3 |0>
        !                         = + a^+_5 a^+_8 a^+_3 |0>
        !                         = - a^+_5 a^+_3 a^+_8 |0>
        !                         = + a^+_3 a^+_5 a^+_8 |0>
        !                         = + | 3 5 8 >

        ! This potential sign change must be taken into account; we do so by
        ! absorbing the sign into the signed population of excips on the
        ! excitor.

        ! Essentially taken from Alex Thom's original implementation.

        ! In:
        !    excitor: bit string of the Slater determinant formed by applying
        !        the excitor to the reference determinant.
        !    excitor_level: excitation level, relative to the determinant f,
        !        of the excitor.  Equal to the number of
        !        annihilation (or indeed creation) operators in the excitor.
        !    f: bit string of determinant to which excitor is
        !       applied to generate a new determinant.
        ! Out:
        !    excitor_sign: sign due to applying the excitor to the
        !       determinant f to form a Slater determinant, i.e. < D_i | a_i D_f >,
        !       which is +1 or -1, where D_i is the determinant formed from
        !       applying the cluster of excitors, a_i, to D_f

        use const, only: i0_end

        integer(i0), intent(in) :: excitor(:)
        integer, intent(in) :: excitor_level
        integer, intent(inout) :: excitor_sign
        integer(i0), intent(in) :: f(:)

        integer(i0) :: excitation(size(excitor))
        integer :: ibasis, ibit, ncreation, nannihilation

        ! Bits involved in the excitation from the reference determinant.
        excitation = ieor(f, excitor)

        nannihilation = excitor_level
        ncreation = excitor_level

        excitor_sign = 1

        ! Obtain sign change by (implicitly) constructing the determinant formed
        ! by applying the excitor to the reference determinant.
        do ibasis = 1, size(excitor)
            do ibit = 0, i0_end
                if (btest(f(ibasis),ibit)) then
                    ! Occupied orbital in reference.
                    if (btest(excitation(ibasis),ibit)) then
                        ! Orbital excited from...annihilate electron.
                        ! This amounts to one fewer operator in the cluster through
                        ! which the other creation operators in the determinant must
                        ! permute.
                        nannihilation = nannihilation - 1
                    else
                        ! Orbital is occupied in the reference and once the
                        ! excitor has been applied.
                        ! Permute the corresponding creation operator through
                        ! the remaining creation and annihilation operators of
                        ! the excitor (which operate on orbitals with a higher
                        ! index than the current orbital).
                        ! If the permutation is odd, then we incur a sign
                        ! change.
                        if (mod(nannihilation+ncreation,2) == 1) &
                            excitor_sign = -excitor_sign
                    end if
                else if (btest(excitation(ibasis),ibit)) then
                    ! Orbital excited to...create electron.
                    ! This amounts to one fewer operator in the cluster through
                    ! which the creation operators in the determinant must
                    ! permute.
                    ! Note due to definition of the excitor, it is guaranteed
                    ! that this is created in the correct place, ie there are no
                    ! other operators in the excitor it needs to be interchanged
                    ! with.
                    ncreation = ncreation - 1
                end if
            end do
        end do

    end subroutine convert_excitor_to_determinant

end module ccmc_data
