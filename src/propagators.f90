module propagators
! For now only includes the wall-Chebyshev propagator
! [TODO] - Maybe move the quasi-Newton propagator here?
use const

implicit none

contains

    subroutine init_chebyshev(qs, sys)
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
        use qmc_data, only: qmc_state_t
        use system, only: sys_t
        use determinant_enumeration, only: enumerate_determinants
        use determinants, only: encode_det
        use proc_pointers, only: sc0_ptr
        use hamiltonian, only: get_hmatel
        use hamiltonian_data, only: hmatel_t

        type(qmc_state_t), intent(inout) :: qs
        type(sys_t), intent(in) :: sys
        
        integer, allocatable :: occ_list_max(:), sym_space_size(:)
        integer(i0), allocatable :: singles_doubles(:,:), f_max(:)
        integer :: ndets, i
        real(p) :: e_max
        type(hmatel_t) :: offdiagel

        allocate(qs%cheby_prop%zeroes(qs%cheby_prop%order))
        allocate(occ_list_max(sys%basis%nbasis))
        call chebyshev_zeroes(qs%cheby_prop, qs%shift(1)) 

        occ_list_max = qs%ref%occ_list0(sys%basis%nbasis:1:-1) ! inverts the occ_list
        call enumerate_determinants(sys, .true., .false., 2, sym_space_size, ndets, singles_doubles, 1, occ_list_max) ! init first
        call enumerate_determinants(sys, .false., .false., 2, sym_space_size, ndets, singles_doubles, 1, occ_list_max) ! store determs
        ! BZ [TODO] - Make sure this works with even selection, which has 1 info_string in tot_string_len
        allocate(f_max(sys%basis%tot_string_len))
        call encode_det(sys%basis, occ_list_max, f_max)

        e_max = sc0_ptr(sys, f_max) + 0.5 ! Arbitrarily shift the upper bound higher to be safe
        do i=1, ndets
            ! BZ [TODO] - Deal with complex systems
            offdiagel = get_hmatel(sys, qs%ref%f0, singles_doubles(:,i))
            e_max = e_max + offdiagel%r
        end do

        qs%cheby_prop%spectral_range(1) = 0.0_p
        qs%cheby_prop%spectral_range(2) = e_max - sc0_ptr(sys, qs%ref%f0)

    end subroutine init_chebyshev

    subroutine update_chebyshev(cheby_prop, shift)

        ! S_i = E_0 + R(1 - cos(i/(m + 1/2)*pi))
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
        end do

    end subroutine update_chebyshev

end module propagators
