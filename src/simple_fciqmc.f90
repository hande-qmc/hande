module simple_fciqmc

! This module contains a very simple, very dumb FCIQMC algorithm.
! This is a serial-only algorithm and uses lots of memory---in particular it
! requires that the Hamiltonian matrix and list of determinants are stored.

! Nonetheless, it is useful for debugging and having a simple algorithm which
! definitely works...

! Based on GHB's ModelFCIQMC code.

use const
use errors

use calc
use determinants
use fciqmc_data

implicit none

contains

    subroutine init_simple_fciqmc()

        use parallel, only: nprocs
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

        ! Allocate main and spawned lists to hold population of walkers.
        allocate(walker_population(ndets), stat=ierr)
        allocate(spawned_walker_population(ndets), stat=ierr)
        ! Zero these.
        walker_population = 0
        spawned_walker_population = 0

        ! Initialise random numbers.
        call random_seed()

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

        write (6,'(1X,a29,'//int_fmt(ref_det)//',/)') 'Reference det is determinant:', ref_det

    end subroutine init_simple_fciqmc

    subroutine do_simple_fciqmc()

        logical :: tVaryShift = .false.
        real(dp) :: shift_damping = 0.050_dp
        integer :: target_particles = 10000
        integer :: ireport, icycle
        integer :: nparticles, nparticles_old

        integer :: iwalker, ipart, j, nspawn, nkill
        real(dp) :: rate
        real(dp) :: r

        do ireport = 1, nreport

            do icycle = 1, ncycles

                ! Zero spawning arrays.
                spawned_walker_population = 0

                ! Consider all walkers.
                spawndie: do iwalker = 1, ndets

                    ! Simulate spawning.
                    ! 1. Consider all particles on current determinant.
                    spawning: do ipart = 1, abs(walker_population(iwalker))
                        ! 2. Simulate spawning by attempting to spawn on all
                        ! connected determinants.
                        do j = 1, ndets

                            ! Can't spawn onto self.
                            if (iwalker == j) cycle
                            ! Can't spawn onto disconnected dets
                            if (hamil(iwalker,j) == 0.0_dp) cycle

                            ! Attempt spawning.
                            rate = abs(Tau*hamil(iwalker,j))
                            nspawn = int(rate)
                            rate = rate - nspawn
                            call random_number(r)
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
                    end do spawning

                    ! Simulate cloning/death.
                    rate = abs(walker_population(iwalker))*tau*(hamil(iwalker,iwalker)-H00-shift)
                    nkill = int(rate)
                    rate = rate - nkill
                    call random_number(r)
                    if (abs(rate) > r) then
                        if (rate > 0.0_dp) then
                            nkill = nkill + 1
                        else
                            nkill = nkill - 1
                        end if
                    end if
                    if (nkill > abs(walker_population(iwalker))) then
                        write (6,*) iwalker, walker_population(iwalker), &
                        abs(walker_population(iwalker))*tau*(hamil(iwalker,iwalker)-H00-shift)
                        call stop_all('do_simple_fciqmc','Trying to create anti-particles.')
                    end if
                    if (walker_population(iwalker) > 0) then
                        walker_population(iwalker) = walker_population(iwalker) - nkill
                    else
                        walker_population(iwalker) = walker_population(iwalker) + nkill
                    end if

                end do spawndie

                ! Annihilation: merge main and spawned lists.
                forall (iwalker = 1:ndets)
                    walker_population(iwalker) = walker_population(iwalker) + spawned_walker_population(iwalker)
                end forall

            end do

            ! Update the shift
            nparticles = sum(abs(walker_population))
            if (tVaryShift) then
                shift = shift - log(real(nparticles,8)/nparticles_old)*shift_damping/(tau*ncycles)
            end if
            nparticles_old = nparticles
            if (nparticles > target_particles) then
                tVaryShift = .true.
            end if
            
            ! Output stats
            write (6,'(i8,f20.10,i8)') ireport*ncycles, shift, nparticles

        end do

    end subroutine do_simple_fciqmc

end module simple_fciqmc
