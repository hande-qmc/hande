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
        real(p) :: spectral_range
        real(p), allocatable :: zeroes(:)

    end type cheb_t

    subroutine init_chebyshev(sys, qs, shift)
        ! Initialises parameters to do with the wall-Chebyshev propagator at the start of the simulation.
        ! In:
        !   sys: system under study.
        !   shift: estimate of the current correlation energy, either the inst. proj. energy or the shift.
        ! In/out:
        !   qs: the qmc_state_t object containing the Chebyshev propagator being initialised.

        ! Initial estimate of the spectral range
        ! Two options: 
        !   1. E_{highest determinant} - E_{HF}
        !   2. Gershgorin circles

        ! Set up initial zeroes
        use determinant_enumeration, only: enumerate_determinants
        use proc_pointers, only: sc0_ptr

        type(cheb_t), intent(inout) :: CHEBPROP
        real(p), intent(in) :: shift
        integer, allocatable :: occ_list_max(:)

        allocate(CHEBPROP%zeroes(CHEBPROP%order))
        allocate(occ_list_max(sys%basis%nbasis))
        call chebyshev_zeroes(CHEBPROP, shift) 
        occ_list_max = qs%ref%occ_list0(sys%basis%nbasis,1,-1) ! inverts the occ_list
        CHEBPROP%spectral_range = sc0_ptr(sys, f_max) - sc0_ptr(sys, reference%f0)


    end subroutine init_wall_chebyshev

    subroutine chebyshev_zeroes(CHEBPROP, shift)

        ! S_i = E_0 + R(1 - cos(i/(m + 1/2)*pi))

        do i = 1, CHEBPROP%order
            CHEBPROP%zeroes(i) = shift + CHEBPROP%spectral_range * (1 - cos(pi*i / (CHEBPROP%order+0.5) ))
        end do

    end subroutine chebyshev_zeroes

end module propagators