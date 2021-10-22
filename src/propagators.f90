module propagators
! For now only includes the wall-Chebyshev propagator
! [TODO] - Maybe move the quasi-Newton propagator here?
use const

implicit none

contains

    type cheb_t
        ! The wall-Chebyshev propagator, where instead of 
        ! (linearly) approximating e^{-\tau H}, we directly approximate \lim_{\tau -> \inf} e^{-\tau H},
        ! which is the "wall function", and we expand the wall function in Chebyshev polynomials.
        ! The Chebyshev polynomials are amenable for use as projectors for the following reasons:
        !   1. They are bound by [-1,1] in the (estimated) spectral range, becoming small near the upper spectral bound
        !   2. They diverge to +\inf as E -> -\inf, meaning the lower spectral bound can be an estimate
        !   3. Sums of up to m-th order Chebyshev polynomials can be written as products of m linear projectors, each involving
        !       the m-th zero of the the original sum.
        !   4. Most importantly they kill off non ground states (arbitrarily, by making m large) faster than the linear projector.
        ! See 10.1021/acs.jctc.6b00639
        integer :: order = 5 ! Default, same as 10.1021/acs.jctc.6b00639
        real(p) :: spectral_range(2) = (/0.0_p, 0.0_p/)
        real(p), allocatable :: zeroes(:)

    end type cheb_t

    subroutine init_chebyshev(sys, qs)
        ! Initialises parameters to do with the wall-Chebyshev propagator at the start of the simulation.
        ! In:
        !   sys: system under study.
        ! In/out:
        !   qs: the qmc_state_t object containing the Chebyshev propagator being initialised.

        ! Initial estimate of the spectral range
        ! Two options: 
        !   1. E_{highest determinant} - E_{HF}
        !   2. Gershgorin circles

        ! Set up initial zeroes
        use determinant_enumeration, only: enumerate_determinants
        use proc_pointers, only: sc0_ptr
        use hamiltonian, only: get_hmatel

        type(qmc_state_t), intent(inout) :: qs
        real(p), intent(in) :: shift
        integer, allocatable :: occ_list_max(:), sym_space_size(:)
        integer(i0), allocatable :: singles_doubles(:,:) ! (tot_string_len,ndets)
        integer :: ndets, i

        allocate(qs%cheby_prop%zeroes(qs%cheby_prop%order))
        allocate(occ_list_max(sys%basis%nbasis))
        call chebyshev_zeroes(qs%cheby_prop, shift) 

        occ_list_max = qs%ref%occ_list0(sys%basis%nbasis,1,-1) ! inverts the occ_list
        call enumerate_determinants(sys, .true., .false., 2, sym_space_size, ndets, singles_doubles, 1, occ_list_max) ! init first
        call enumerate_determinants(sys, .false., .false., 2, sym_space_size, ndets, singles_doubles, 1, occ_list_max) ! store determs
        ! BZ [TODO] - Make sure this works with even selection, which has 1 info_string in tot_string_len

        e_max = sc0_ptr(sys, f_max) + 0.5 ! Arbitrarily shift the upper bound higher to be safe
        do i=1, ndets
            e_max += get_hmatel(sys, qs%ref%f0, singles_doubles[i])
        end do

        qs%cheby_prop%spectral_range(2) = e_max - sc0_ptr(sys, reference%f0)

    end subroutine init_wall_chebyshev

    subroutine update_chebyshev(cheby_prop, shift)

        ! S_i = E_0 + R(1 - cos(i/(m + 1/2)*pi))
        ! In:
        !   shift: the current shift, needed to update the zeroes.
        ! In/out:
        !   cheby_prop: the cheb_t object containing the Chebyshev propagator.
        real(p), intent(in) :: shift
        type(cheb_t), intent(inout) :: cheby_prop
        integer :: i

        cheby_prop%spectral_range(1) = shift

        do i = 1, cheby_prop%order
            cheby_prop%zeroes(i) = shift + (cheby_prop%spectral_range(2)-cheby_prop%spectral_range(1))/2 &
                                   * (1 - cos(pi*i / (cheby_prop%order+0.5) )) ! BZ [TODO] - explain the factor of 2
        end do

    end subroutine update_chebyshev

end module propagators