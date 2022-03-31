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

    pure function particle_number(sys, beta, mu, ispin) result (nav)

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
        !    mu: initial guess for the chemical potential.
        !    ispin: integer to account for odd/even ordering of alpha/beta
        !       spin orbitals. Set to 1 for alpha spins, 2 for beta spins.
        ! Returns:
        !    nav: expected number of particles.

        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        real(p), intent(in) :: beta, mu
        integer, intent(in) :: ispin

        real(p) :: nav

        integer :: iorb
        real(p) :: ff, ek

        nav = 0.0_p

        do iorb = ispin, sys%basis%nbasis, 2
            ek = sys%basis%basis_fns(iorb)%sp_eigv
            ff = fermi_factor(ek, mu, beta)
            nav = nav + ff
        end do

    end function particle_number

    pure function deriv_particle_number(sys, beta, mu, ispin) result(nav_deriv)

        ! Derivative of the expected number of particles with respect to the
        ! chemical potential, i.e.,
        ! Evaluate:

        ! d/dmu N |_{mu=mu'} = d/dmu sum_k f_k
        !                    = beta sum_k e^{beta(ek-mu')}/(e^{beta(ek-mu')}+1)^2
        !                    = beta \sum_k 1 / (2(cosh(beta(ek-mu'))+1))
        ! where we go from the second to the third line for (numerical) stability purposes.

        ! In:
        !    basis: type containing information on the single-particle basis
        !        set.
        !    beta: inverse temperature.
        !    mu: chemical potential.
        !    ispin: integer to account for odd/even ordering of alpha/beta
        !       spin orbitals. Set to 1 for alpha spins, 2 for beta spins.
        ! Returns:
        !    nav_deriv: derivative of the N wrt mu.

        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        real(p), intent(in) :: beta, mu
        real(p) :: nav_deriv
        integer, intent(in) :: ispin

        integer :: iorb
        real(p) :: ek

        nav_deriv = 0.0_p

        do iorb = ispin, sys%basis%nbasis, 2
            ek = sys%basis%basis_fns(iorb)%sp_eigv
            nav_deriv = nav_deriv + beta/(2*(cosh(beta*(ek-mu))+1))
        end do

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

        use system, only: sys_t, ueg
        use errors, only: stop_all
        use parallel

        type(sys_t), intent(in) :: sys
        real(p), intent(in) :: beta

        integer :: it, max_it, ispin, nels(1:2)
        real(p) :: nav, nav_deriv, scale_factor, mu(1:2), mu_old(1:2)
        logical :: conv(1:2)

        max_it = 10000 ! In case it doesn't converge.
        it = 0

        if (sys%nvirt_alpha == 0 .or. sys%nvirt_beta == 0) then
            call stop_all('find_chem_pot', 'Completely filled spin channels are not &
                                            supported, please increase the basis set.')
        end if

        ! The alpha and beta channels have their chemical
        ! potentials stored seperatly in the mu array.
        ! In the basis_fns array:
        ! alpha orbitals start at index 1 and appear every 2
        ! beta orbitals start at index 2 and appear every 2
        ! Hence the use of ispin of 1 and 2.
        nels = (/ sys%nalpha, sys%nbeta /)
        do ispin = 1, 2
            if (nels(ispin) == 0) then
                ! If we have no particles we should have the vacuum chemical
                ! potential which is zero.
                mu_old(ispin) = 0.0_p
                conv(ispin) = .true.
            else
                ! Otherwise we want to start at the T = 0 fermi energy
                ! which is typically 0.5*(HOMO + LUMO)
                ! Note that we took care of the case of no HOMO and no LUMO above
                mu_old(ispin) = sys%basis%basis_fns(2*nels(ispin)-(2-ispin))%sp_eigv
                mu_old(ispin) = (mu_old(ispin) + sys%basis%basis_fns(2*nels(ispin)+ispin)%sp_eigv)/2.0_p
                conv(ispin) = .false.
            end if
        end do

        if (sys%system == ueg) then
            ! At high density and large basis sets (at higher temperatures) the chemical
            ! potential can grow to be quite large, so rescale the threshold by an appropriate
            ! amount to avoid numerical precision issues.
            scale_factor = sys%ueg%ef
        else
            ! [todo] - Update if necessary.
            scale_factor = 1.0_p
        end if

        do while (it < max_it)
            do ispin = 1, 2
                if (.not. conv(ispin)) then
                    nav = particle_number(sys, beta, mu_old(ispin), ispin)
                    nav_deriv = deriv_particle_number(sys, beta, mu_old(ispin), ispin)
                    mu(ispin) = mu_old(ispin) - (nav-nels(ispin))/nav_deriv
                    if (abs(mu(ispin)-mu_old(ispin)) / scale_factor < depsilon) then
                        ! Found the root!
                        conv(ispin) = .true.
                    else
                        mu_old(ispin) = mu(ispin)
                    end if
                end if
            end do
            it = it + 1
            if (conv(1) .and. conv(2)) exit
        end do

        if (it == max_it) then
            call stop_all('find_chem_pot', 'Root not found, check interval.')
        end if

    end function find_chem_pot

end module chem_pot
