module simple_fciqmc

! This module contains a very simple, very dumb, self-contained FCIQMC algorithm.
! This is a serial-only algorithm and uses lots of memory---in particular it
! requires that the Hamiltonian matrix and list of determinants are stored.

! Nonetheless, it is useful for debugging and having a simple algorithm which
! definitely works...

! Based on GHB's ModelFCIQMC code.

use const
use dSFMT_interface
use errors

use calc
use determinants
use fciqmc_data

implicit none

contains

    subroutine init_simple_fciqmc()

        use parallel, only: nprocs, parent
        use utils, only: int_fmt

        use diagonalisation, only: generate_hamil

        integer :: ierr
        integer :: i, j, ref_det

        if (nprocs > 1) call stop_all('init_simple_fciqmc','Not a parallel algorithm.')

        ! Find and set information about the space.
        call set_spin_polarisation(ms_in)
        call find_sym_space_size()

        ! Find all determinants with desired spin and symmetry.
        call enumerate_determinants(sym_in)

        ! Set up hamiltonian matrix.
        call generate_hamil(distribute_off)
        ! generate_hamil fills in only the lower triangle.
        ! fill in upper triangle for easy access.
        do i = 1,ndets
            do j = i+1, ndets
                hamil(j,i) = hamil(i,j)
            end do
        end do

        write (6,'(1X,a13,/,1X,13("-"),/)') 'Simple FCIQMC'
        write (6,'(1X,a53,1X)') 'Using a simple (but correct) serial FCIQMC algorithm.'
        write (6,'(1X,a137)') 'Enumeration of the determinant list and evaluation of &
                              &the Hamiltonian matrix for the given symmetry block and &
                              &spin polarization required.'
        write (6,'(1X,a104,/)') 'This is slow and memory demanding: consider using the &
                                &fciqmc option instead of the simple_fciqmc option.'
        write (6,'(1X,a46,'//int_fmt(sym_in,1)//',1X,a9,'//int_fmt(ms_in,1)//',a1,/)') &
            'Considering determinants belonging to symmetry',sym_in,'with spin',ms_in,"."

        ! Allocate main and spawned lists to hold population of walkers.
        allocate(walker_population(ndets), stat=ierr)
        allocate(spawned_walker_population(ndets), stat=ierr)
        ! Zero these.
        walker_population = 0
        spawned_walker_population = 0

        ! Initialise random numbers.
        call dSFMT_init(5234)

        ! Now we need to set the reference determinant.
        ! We choose the determinant with the lowest Hamiltonian matrix element.
        ref_det = 1
        do i = 2, ndets
            if (hamil(i,i) < hamil(ref_det, ref_det)) then
                ref_det = i
            end if
        end do

        H00 = hamil(ref_det,ref_det)
        walker_population(ref_det) = 1

        write (6,'(1X,a29,1X)',advance='no') 'Reference determinant, |D0> ='
        call write_det(dets_list(:,ref_det), new_line=.true.)
        write (6,'(1X,a16,f20.12)') 'E0 = <D0|H|D0> =',H00
        write (6,'(/,1X,a68,/)') 'Note that FCIQMC calculates the correlation energy relative to |D0>.'

    end subroutine init_simple_fciqmc

    subroutine do_simple_fciqmc()

        integer :: ireport, icycle, iwalker, ipart
        integer :: nparticles, nparticles_old

        write (6,'(1X,a12,7X,a13,2X,a11)') '# iterations','Instant shift','# particles'

        do ireport = 1, nreport

            do icycle = 1, ncycles

                ! Zero spawning arrays.
                spawned_walker_population = 0

                ! Consider all walkers.
                do iwalker = 1, ndets

                    ! Simulate spawning.
                    do ipart = 1, abs(walker_population(iwalker))
                        ! Attempt to spawn from the current particle onto all
                        ! connected determinants.
                        call attempt_spawn(iwalker)
                    end do

                    call simple_death(iwalker)

                end do

                call simple_annihilation()

            end do

            ! Update the shift
            nparticles = sum(abs(walker_population))
            if (vary_shift) then
                call update_shift(nparticles_old, nparticles, ncycles)
            end if
            nparticles_old = nparticles
            if (nparticles > target_particles) then
                vary_shift = .true.
            end if
            
            ! Output stats
            write (6,'(5X,i8,f20.10,2X,i11)') ireport*ncycles, shift, nparticles

        end do

        write (6,'(/,1X,a13,f22.12)') 'final_shift =',shift
        write (6,'(1X,a12,1X,f22.12)') 'E0 + shift =',shift+H00

    end subroutine do_simple_fciqmc

    subroutine attempt_spawn(iwalker)

        ! Simulate spawning part of FCIQMC algorithm.
        ! In:
        !    iwalker: walker whose particles attempt to clone/die.

        integer, intent(in) :: iwalker

        integer :: j, nspawn
        real(dp) :: rate
        real(dp) :: r

        ! Simulate spawning by attempting to spawn on all
        ! connected determinants.
        do j = 1, ndets

            ! Can't spawn onto self.
            if (iwalker == j) cycle
            ! Can't spawn onto disconnected dets
            if (hamil(iwalker,j) == 0.0_dp) cycle

            ! Attempt spawning.
            ! Spawn with probability tau|K_ij|.
            ! As K_ij = H_ij for off-diagonal elements, we can just use the
            ! stored Hamiltonian matrix directly.
            rate = abs(Tau*hamil(iwalker,j))
            nspawn = int(rate)
            rate = rate - nspawn
            r = genrand_real2()
            if (rate > r) nspawn = nspawn + 1

            ! Create particles.
            if (hamil(iwalker,j) > 0.0_dp) then
                ! Flip child sign.
                if (walker_population(iwalker) < 0) then
                    ! Positive offspring.
                    spawned_walker_population(j) = spawned_walker_population(j) + nspawn
                else
                    spawned_walker_population(j) = spawned_walker_population(j) - nspawn
                end if
            else
                ! Same sign as parent.
                if (walker_population(iwalker) > 0) then
                    ! Positive offspring.
                    spawned_walker_population(j) = spawned_walker_population(j) + nspawn
                else
                    spawned_walker_population(j) = spawned_walker_population(j) - nspawn
                end if
            end if

        end do

    end subroutine attempt_spawn

    subroutine simple_death(iwalker)

        ! Simulate cloning/death part of FCIQMC algorithm.
        ! In:
        !    iwalker: walker whose particles attempt to clone/die.

        integer, intent(in) :: iwalker

        integer :: nkill
        real(dp) :: rate
        real(dp) :: r

        ! A particle dies with probability, p_d, given by
        !  p_d = tau(K_ii _ S)
        ! where tau is the timestep, S is the shift and K_ii is
        !  K_ii =  < D_i | H | D_i > - E_0
        ! We store the Hamiltonian matrix rather than the K matrix.
        ! It is efficient to allow all particles on a given determinant to
        ! attempt to die in one go (like lemmings) in a stochastic process.
        rate = abs(walker_population(iwalker))*tau*(hamil(iwalker,iwalker)-H00-shift)
        ! Number to definitely kill.
        nkill = int(rate)
        rate = rate - nkill

        ! Additional stochasitic death?
        r = genrand_real2()
        if (abs(rate) > r) then
            if (rate > 0.0_dp) then
                nkill = nkill + 1
            else
                nkill = nkill - 1
            end if
        end if

        ! Don't allow creation of anti-particles in simple_fciqmc.
        if (nkill > abs(walker_population(iwalker))) then
            write (6,*) iwalker, walker_population(iwalker), &
            abs(walker_population(iwalker))*tau*(hamil(iwalker,iwalker)-H00-shift)
            call stop_all('do_simple_fciqmc','Trying to create anti-particles.')
        end if

        ! Update walker populations.
        ! Particle death takes the walker_population closer to 0...
        ! (and similarly if cloning (ie nkill is negative) then the
        ! walker_population should move away from 0...)
        if (walker_population(iwalker) > 0) then
            walker_population(iwalker) = walker_population(iwalker) - nkill
        else
            walker_population(iwalker) = walker_population(iwalker) + nkill
        end if

    end subroutine simple_death

    subroutine simple_annihilation()

        ! Annihilation: merge main and spawned lists.

        ! This is especially easy as we store the walker populations for all
        ! determinants for both the main and spawned lists so it just amounts to
        ! adding the two arrays together,

        walker_population = walker_population + spawned_walker_population

    end subroutine simple_annihilation

    subroutine update_shift(nparticles_old, nparticles,nupdate_steps)

        ! Update the shift according to:
        !  shift(beta) = shift(beta-A*tau) - xi*log(N_w(tau)/N_w(beta-A*tau))/(A*tau)
        ! where
        !  * shift(beta) is the shift at imaginary time beta;
        !  * A*tau is the amount of imaginary time between shift-updates (=# of
        !    Monte Carlo cycles between updating the shift);
        !  * xi is a damping factor (0.05-0.10 is appropriate) to prevent large fluctations;
        !  * N_w(beta) is the total number of particles at imaginary time beta.
        ! In:
        !    nparticles_old: N_w(beta-A*tau).
        !    nparticles: N_w(beta).
        !    

        integer, intent(in) :: nparticles_old, nparticles, nupdate_steps

        ! This should be changed into an input option when necessary.
        real(dp) :: shift_damping = 0.050_dp

        shift = shift - log(real(nparticles,8)/nparticles_old)*shift_damping/(tau*nupdate_steps)

    end subroutine update_shift

end module simple_fciqmc
