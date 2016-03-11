module ct_fciqmc

! Evolve the walker population using a continuous time algorithm (i.e. jump
! directly to the next spawning event without a timestep).

use fciqmc_data
use const, only: p, int_64

implicit none

contains


    subroutine do_ct_fciqmc(sys, qmc_in, restart_in, reference, load_bal_in, annihilation_flags, matel)

        ! In:
        !    sys: system being studied
        !    qmc_in: Input options relating to QMC methods.
        !    restart_in: input options for HDF5 restart files.
        !    reference: reference determinant.
        !    matel: off-diagonal Hamiltonian matrix element (ignoring sign due
        !       to permutations).  Either U (Bloch orbitals) or
        !       t (atomic/real-space orbitals).
        !    load_bal_in: input options for load balancing.
        !    annihilation_flags: calculation specific annihilation flags.

        use system, only: sys_t
        use interact

        use qmc_data, only: qmc_in_t, restart_in_t, reference_t, load_bal_in_t, annihilation_flags_t
        use errors, only: stop_all

        type(sys_t), intent(in) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        type(restart_in_t), intent(in) :: restart_in
        type(reference_t), intent(in) :: reference
        real(p), intent(in) :: matel ! either U or t, depending whether we are working in the real or k-space
        type(load_bal_in_t), intent(in) :: load_bal_in
        type(annihilation_flags_t), intent(in) :: annihilation_flags

        call stop_all('do_ct_fciqmc', 'Awaiting rewrite from JSS.')

    end subroutine do_ct_fciqmc

    function timestep(rng, R) result(dt)

        ! In:
        !    R: \sum_i < D | H - E_0 - S | D_i >, the sum of all non-zero matrix
        !       elements connected to the current determinant, D.
        ! In/Out:
        !    rng: random number generator.
        ! Returns:
        !    dt: the (stochastic) time which elapses before the next spawning
        !        event.

        use dSFMT_interface, only: dSFMT_t, get_rand_close_open

        real(p) :: dt
        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(in)  :: R

        dt = -(1.0_p/R)*log(real(get_rand_close_open(rng),p))

    end function timestep


end module ct_fciqmc
