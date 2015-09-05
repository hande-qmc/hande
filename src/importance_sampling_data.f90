module importance_sampling_data

! Data for performing importance sampling.
! Currently this is rather nascent but, should it prove useful, be easily extended
! to other trial wavefunctions, systems and importance sampling methods.

use const, only: p

implicit none

! Types of wavefunctions implemented:

! For the Heisenberg model, several different trial functions can be used in the
! projected energy estimator. Only a single determinant can be used for other systems.
enum, bind(c)
    enumerator :: single_basis = 0
    enumerator :: neel_singlet
end enum

! Types of Hamiltonian transformations implemented:

! For the Heisenberg model, a guiding function may be used,
! |psi_G> = \sum_{i} a_i |psi_i>, so that the new Hamiltonian matrix elements are
! H_ij^new = (a_i*H_ij)/a_j. This is just importance sampling.
! Possible |psi_G>:
! NOTE: DMQMC importance sampling is handled in dmqmc_data.
enum, bind(c)
    ! a_i = 1.
    enumerator :: no_guiding
    ! Note that when we use the Neel singlet state as a guiding function, it must also
    ! be used as the trial function in calculating the projected energy.
    enumerator :: neel_singlet_guiding
end enum

! === WARNING ===
! trial_t%wfn_dat has variable bounds depending upon the wavefunction in use.  In
! particular, if passing trial_t%wfn_dat into a procedure, care should be taken
! to ensure the lower bounds are correctly defined for the desired wavefunction.

! Trial wavefunction information.
type trial_t
    ! Trial wavefunction in use.  Must be a member of the trial function enum.
    integer :: wfn = single_basis
    ! Guiding wavefunction/Hamiltonian transformation in use.  Note: not all combinations
    ! of wfn and guide are implemented!
    integer :: guide = no_guiding
    ! Wavefunction-specific data.
    ! single_basis:
    !     Not used.
    ! neel_singlet:
    !     Neel singlet amplitudes.  wfn_dat(n) is the amplitude
    !     of a state with n up spins in the first sublattice.
    !     dimensions: (-1:nspins/2+1) and wfn_dat(-1) = wfn_dat(nspins/2+1) = 0.
    real(p), allocatable :: wfn_dat(:)
end type trial_t

end module importance_sampling_data
