module chem_pot

! Module to evaluate the chemical potential at a given temperature from an underlying single-particle basis set.

use const, only: p, depsilon

implicit none

contains

    pure function fermi_factor(ek, mu, beta) result(f_k)

        ! Calculate the Fermi factor, f_k = 1/(e^{beta(e_k-mu)}+1).

        ! In:
        !    ek: single particle eigen value.
        !    mu: chemical potential.
        !    beta: inverse temperature.
        ! Returns:
        !    fk: Fermi factor for single-particle state ek.

        real(p), intent(in) :: ek, mu, beta

        real(p) :: f_k

        f_k = 1.0_p / (exp(beta*(ek-mu))+1)

    end function fermi_factor

    pure function particle_number(sys, beta, mu) result (nav)

        ! Evaluate the expected number of particles for given chemical potential
        ! mu at temperature T = 1 / beta.

        ! In the grand canonical ensemble the expected number of particles is
        ! given as:
        !
        !     N = \sum_k f_k(beta, mu),
        !
        ! for f_k a fermi factor defined above.

        ! In:
        !    basis: type containing information on the single-particle basis
        !        set.
        !    beta: inverse temperature.
        !    mu_init: initial guess for the chemical potential.
        ! Returns:
        !    nav: expected number of particles.

        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        real(p), intent(in) :: beta, mu

        real(p) :: nav

        integer :: iorb
        real(p) :: ff, ek

        nav = 0.0_p

        do iorb = 1, sys%basis%nbasis, 2
            ek = sys%basis%basis_fns(iorb)%sp_eigv
            ff = fermi_factor(ek, mu, beta)
            nav = nav + ff
        end do

        ! Unpolarised?
        if (sys%ms == 0) nav = 2.0_p * nav

    end function particle_number

    pure function deriv_particle_number(sys, beta, mu) result(nav_deriv)

        ! Derivative of the expected number of particles with respect to the
        ! chemical potential, i.e.,
        ! Evaluate:

        ! d/dmu N |_{mu=mu'} = d/dmu sum_k f_k
        !                    = beta sum_k e^{beta(ek-mu')}/(e^{beta(ek-mu')}+1)^2

        ! In:
        !    basis: type containing information on the single-particle basis
        !        set.
        !    beta: inverse temperature.
        !    mu: chemical potential.
        ! Returns:
        !    nav_deriv: derivative of the N wrt mu.

        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        real(p), intent(in) :: beta, mu
        real(p) :: nav_deriv

        integer :: iorb
        real(p) :: ff, ek

        nav_deriv = 0.0_p

        do iorb = 1, sys%basis%nbasis, 2
            ek = sys%basis%basis_fns(iorb)%sp_eigv
            ff = fermi_factor(ek, mu, beta)
            nav_deriv = nav_deriv + beta*exp(beta*(ek-mu))*ff**2.0_p
        end do

        ! Unpolarised?
        if (sys%ms == 0) nav_deriv = 2.0_p * nav_deriv

    end function deriv_particle_number

    function find_chem_pot(sys, beta) result (mu)

        ! Find the chemical potential (mu) at temperature T = 1/beta.

        ! In the grand canonical ensemble the expected number of particles is
        ! given as:
        !
        !     N = \sum_k f_k(beta, mu),
        !
        ! for f_k a fermi factor defined above. To find the mu which would
        ! give us the desired N we need to solve for:
        !
        !   N - \sum_k f_k(beta,mu) = 0.
        !
        ! This is done using the Newton-Raphson method (see, for example,
        ! https://en.wikipedia.org/wiki/Newton%27s_method.)

        ! In:
        !    sys: system under consideration.
        !    beta: inverse temperature.
        ! Returns:
        !    mu: chemical potential of system at inverse temperature beta.

        use system, only: sys_t
        use errors, only: stop_all

        type(sys_t), intent(in) :: sys
        real(p), intent(in) :: beta

        integer :: it, max_it
        real(p) :: mu, mu_old, nav, nav_deriv

        max_it = 10000 ! In case it doesn't converge.

        it = 0

        ! At zero temperature the chemical potential should be close to the
        ! energy of the highest occupied orbital (Fermi Energy).
        ! This is typically a good guess for root finding as the chemical
        ! potential is a monotonically decreasing function of temperature, and
        ! its value at T=0 is the Fermi energy (in the TDL).
        ! The following takes care of the spin polarised and unpolarised case given
        ! how we enumerate basis functions. We assume we only deal with completely
        ! polarised or unpolarised systems
        mu_old = sys%basis%basis_fns(2*sys%nalpha)%sp_eigv

        if (sys%ms /= sys%nel .and. sys%ms /= 0) then
            call stop_all('find_chem_pot', 'Spin polarisation requested not supported - please implement.')
        end if

        do while (it < max_it)
            nav = particle_number(sys, beta, mu_old)
            nav_deriv = deriv_particle_number(sys, beta, mu_old)
            mu = mu_old - (nav-sys%nel)/nav_deriv
            if (abs(mu-mu_old) < depsilon) then
                ! Found the root!
                exit
            else
                it = it + 1
                mu_old = mu
            end if
        end do

        if (it == max_it) then
            call stop_all('find_chem_pot', 'Root not found, check interval.')
        end if

    end function find_chem_pot

end module chem_pot
