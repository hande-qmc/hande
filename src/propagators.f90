module propagators
! For now only includes the wall-Chebyshev propagator
! [TODO] - Maybe move the quasi-Newton propagator here?
use const

implicit none

contains

    subroutine init_chebyshev(sys, qmc_in, qs)
        ! Initialises parameters to do with the wall-Chebyshev propagator at the start of the simulation.
        ! This sets up the initial estimate of the spectral range of the Hamiltonian by Gershgorin circle theorem,
        ! and consequently sets up the zeroes of the polynomial expansion.
        ! In:
        !   sys: system under study.
        ! In/out:
        !   qs: the qmc_state_t object containing the Chebyshev propagator being initialised.
        use qmc_data, only: qmc_state_t, qmc_in_t
        use system, only: sys_t
        use determinant_enumeration, only: enumerate_determinants
        use determinants, only: encode_det
        use proc_pointers, only: sc0_ptr
        use hamiltonian, only: get_hmatel
        use hamiltonian_data, only: hmatel_t
        use parallel, only: parent

        type(sys_t), intent(in) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        type(qmc_state_t), intent(inout) :: qs
        
        integer, allocatable :: occ_list_max(:), sym_space_size(:)
        integer(i0), allocatable :: singles_doubles(:,:), f_max(:)
        integer :: ndets, i
        real(p) :: e_max
        type(hmatel_t) :: offdiagel
        integer :: iunit
        real :: t1, t2

        iunit = 6

        qs%cheby_prop%using_chebyshev = qmc_in%chebyshev
        qs%cheby_prop%icheb = 1
        ! Default is 1
        qs%cheby_prop%order = qmc_in%chebyshev_order
        
        if (qs%cheby_prop%using_chebyshev) then
            call cpu_time(t1)

            allocate(qs%cheby_prop%zeroes(qs%cheby_prop%order))
            allocate(qs%cheby_prop%weights(qs%cheby_prop%order))
            allocate(occ_list_max(sys%nel))
            ! the highest occ_list is [nbasis-nel:nbasis]
            do i=1, sys%nel
                occ_list_max(i) = sys%basis%nbasis - sys%nel + i
            end do
            call enumerate_determinants(sys, .true., .false., 2, sym_space_size, ndets, singles_doubles,&
                                        sys%symmetry, occ_list_max) ! init first
            call enumerate_determinants(sys, .false., .false., 2, sym_space_size, ndets, singles_doubles,&
                                        sys%symmetry, occ_list_max) ! store determs
            ! BZ [TODO] - Make sure this works with even selection, which has 1 info_string in tot_string_len
            allocate(f_max(sys%basis%tot_string_len))
            call encode_det(sys%basis, occ_list_max, f_max)

            e_max = sc0_ptr(sys, f_max) + 0.5 ! Arbitrarily shift the upper bound higher to be safe
            do i=1, ndets
                ! BZ [TODO] - Deal with complex systems
                offdiagel = get_hmatel(sys, f_max, singles_doubles(:,i))
                e_max = e_max + offdiagel%r
            end do

            qs%cheby_prop%spectral_range(1) = 0.0_p
            qs%cheby_prop%spectral_range(2) = e_max - qs%ref%H00
            call update_chebyshev(qs%cheby_prop, 0.0_p)

            call cpu_time(t2)
            if (parent) then
                write(iunit, '(1X, "# Starting the wall-Chebyshev propagator initialisation.")')
                write(iunit, '(1X, "Chebyshev order:", 1X, I0)') qs%cheby_prop%order
                write(iunit, '(1X, "Initial estimate of spectral range:", 1X, "[", I0, ",", ES15.8, "]")') &
                    0, qs%cheby_prop%spectral_range(2)
                write(iunit, '(1X, "Initial zeroes of the Chebyshev polynomial and weights:")')
                write(iunit, '(1X, "i", 6X, "S_i", 10X, "1/(S_i-E_0)")')
                do i = 1, qs%cheby_prop%order
                    write(iunit, '(1X, I0, 1X, ES15.8, 1X, ES15.8)') &
                        i, qs%cheby_prop%zeroes(i), qs%cheby_prop%weights(i)
                end do
                write(iunit, '(1X, "# Finishing the wall-Chebyshev propgator initialisation, time taken:", 1X, ES15.8)') &
                    t2-t1
            end if
        else
            allocate(qs%cheby_prop%weights(1))
            qs%cheby_prop%weights(1) = 1.0_p ! Similar to QN, the hmatel is always scaled to save on too many checks (branching)
        end if

    end subroutine init_chebyshev

    subroutine update_chebyshev(cheby_prop, shift)
        ! The zeroes, S_i, of the m-th order Chebyshev expansion of the wall function are given by
        ! S_i = E_0 + R(1 - cos(i/(m + 1/2)*pi))
        ! Where E_0 is the shift (estimate of the lowest eigenvalue).
        ! In:
        !   shift: the current shift, needed to update the zeroes.
        ! In/out:
        !   cheby_prop: the cheb_t object containing the Chebyshev propagator.
        use qmc_data, only: cheb_t
        
        real(p), intent(in) :: shift
        type(cheb_t), intent(inout) :: cheby_prop
        integer :: i

        cheby_prop%spectral_range(1) = shift

        do i = 1, cheby_prop%order
            cheby_prop%zeroes(i) = shift + (cheby_prop%spectral_range(2)-cheby_prop%spectral_range(1))/2 &
                                   * (1 - cos(pi*i / (cheby_prop%order+0.5) )) ! BZ [TODO] - explain the factor of 2
            ! This takes care of the absence of the minus sign in front of <D_m|g_wall_ch|D_0>
            cheby_prop%weights(i) = 1/(cheby_prop%zeroes(i)-shift)
        end do

    end subroutine update_chebyshev

end module propagators
