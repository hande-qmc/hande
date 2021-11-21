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

        ! Here we briefly document the wall-Chebyshev propagator:
        ! The true projector that turns any trial wavefunction that are nonorthogonal to the ground state wavefunction has action
        ! g_GS|\psi_trial> ~ |\Psi_CC>
        ! and has the form of 
        ! g_GS = \delta(H-E_CC)
        ! This is the wall function, and the infinite time limit of the exponential projector:
        ! lim_{t->\infty} e^(-\tau H)
        ! Knowing the wall function exactly is equivalent to solving the ground state of H, but we can approximate it with a
        ! suitable polynomial expansion that has the following properties:
        !   1. It is bound by [-1,1] in the (estimated) spectral range, becoming small near the upper spectral bound
        !   2. It diverges to +\infty as E -> -\infty, meaning the lower spectral bound can be an estimate
        ! The Chebyshev polynomials satisfy these conditions, with an additional, very fortunate property
        ! that the sum of all Chebyshev polynomials from orders 1 through m can be written as a product of m linear functions.
        ! This means that we can reuse our machinery of linear propagators to exactly reproduce the wall-Chebyshev projection.
        ! To be exact, the m-th order Chebyshev expansion of the wall function can be written as
        ! g_wall-Ch^m = \prod_{i=1}^m \frac{H-a_i}{S-a_i}
        !   where m is a_i's are the zeroes of the polynomial, and S is the arbitrary shift / estimate of the lower bound.
        ! The zeroes are given in a closed form:
        !   a_i = S + R(1 - cos(i/(m + 1/2)*pi))
        !   where R is the spectral radius (E_{N-1} - E_0), E_{N-1} can be estimated with the Gershgorin disc theorem (see below)    
        !   and E_0 is estimated by S
        ! The action of the g_wall-Ch(^m) is:
        !                         g_wall-Ch|\Psi^(n,0)> = |\Psi^(n+1,0)>
        ! \prod_{i=1}^m \frac{H-a_i}{S-a_i}|\Psi^(n,0)> = |\Psi^(n+1,0)>
        ! and this can be written iteratively as
        !               \frac{H-a_1}{S-a_1}|\Psi^(n,0)> = |\Psi^(n,1)>
        !               \frac{H-a_2}{S-a_2}|\Psi^(n,1)> = |\Psi^(n,2)>
        !               and so on...
        ! We project the equations in classic CC fashion:
        !         <D_m|\Psi> = <D_m|\frac{H-a_i}{S-a_i}|\Psi> 
        !              t~_m^ = -\frac{1}{a_i-S} (\sum_{n \neq m} H_mn t~_n + (H_mm - a_i)t~_m)
        !   t~_m - t_m + t_m = -\frac{1}{a_i-S} (\sum_{n \neq m} H_mn t~_n + (H_mm - a_i)t~_m)
        !                t_m = t_m -\frac{1}{a_i-S} (\sum_{n \neq m} H_mn t~_n + (H_mm - S)t~_m)
        ! the last equation strongly resembles the original update equation
        !                t_m = t_m - dt (\sum_{n \neq m) H_mn t~_n + (H_mm - S)t~_m)
        ! but just with all amplitudes scaled by \frac{1}{a_i-S} (which we term "Chebyshev weights"), and dt (tau) set to unity.

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
            call highest_det(sys, occ_list_max)
            call enumerate_determinants(sys, .true., .false., 2, sym_space_size, ndets, singles_doubles,&
                                        sys%symmetry, occ_list_max) ! init first
            call enumerate_determinants(sys, .false., .false., 2, sym_space_size, ndets, singles_doubles,&
                                        sys%symmetry, occ_list_max) ! store determs
            ! BZ [TODO] - Make sure this works with even selection, which has 1 info_string in tot_string_len
            allocate(f_max(sys%basis%tot_string_len))
            call encode_det(sys%basis, occ_list_max, f_max)

            e_max = 0.5 ! Arbitrarily shift the upper bound higher to be safe
            do i=1, ndets
                ! BZ [TODO] - Deal with complex systems
                ! Note 'offdiagel' here actually includes the diagonal element, as enumerate_determinant returns it
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

    subroutine highest_det(sys, occ_list_max)
        ! This returns the occupation list of the highest determinant of the same spin symmetry. Getting the spatial symmetry 
        ! right too would be nice but wouldn't be necessary for our estimate: we just need an upper bound.
        ! BZ [TODO] perhaps this belongs in one of the determinant source files
        ! This currently only supports read_in systems (extension to other systems should be trivial)
        ! In:
        !   sys: the system under study. Contains basis set information.
        ! Out:
        !   occ_list_max: the occupation list of the highest determinant with the same Ms.
        !                 e,g. for the H4 molecule in STO-3G this would be (/5,6,7,8/)
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer, allocatable, intent(out) :: occ_list_max(:)
        integer :: i_alpha, i_beta
        
        allocate(occ_list_max(sys%nel))

        do i_alpha=1, sys%nalpha
            ! We assume that # of alpha basis fn is equal to # of beta ones
            occ_list_max(i_alpha) = sys%basis%nbasis - (i_alpha*2 - 1) ! (/nbasis-1, nbasis-3, nbasis-5, .../)
        end do

        do i_beta=1, sys%nbeta
            occ_list_max(sys%nalpha+i_beta) = sys%basis%nbasis - (i_beta-1)*2 ! (/nbasis, nbasis-2, .../)
        end do
        ! Note that occ_list_max is not sorted, as this is eventually passed into encode_det where sorting is not necessary
    end subroutine highest_det
            

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
