module simple_fciqmc

! This module contains a very simple, very dumb, self-contained FCIQMC algorithm.
! This is a serial-only algorithm and uses lots of memory---in particular it
! requires that the Hamiltonian matrix and list of determinants are stored.

! Nonetheless, it is useful for debugging and having a simple algorithm which
! definitely works...

use const
use dSFMT_interface
use errors

use calc
use fciqmc_data

implicit none

contains

    subroutine init_simple_fciqmc()

        ! Initialisation for the simple fciqmc algorithm.
        ! Setup the list of determinants in the space, calculate the relevant
        ! symmetry block of the Hamiltonian matrix, initialise the RNG, allocate
        ! the required memory for the list of walkers and set the initial
        ! walker.

        use parallel, only: nprocs, parent
        use checking, only: check_allocate
        use utils, only: int_fmt

        use determinant_enumeration
        use diagonalisation, only: generate_hamil
        use fciqmc_restart, only: read_restart

        integer :: ierr
        integer :: i, j

        if (nprocs > 1) call stop_all('init_simple_fciqmc','Not a parallel algorithm.')

        ! Find and set information about the space.
        call set_spin_polarisation(ms_in)
        if (allocated(occ_list0)) then
            call enumerate_determinants(.true., .false., occ_list0=occ_list0)
        else
            call enumerate_determinants(.true., .false.)
        end if

        ! Find all determinants with desired spin and symmetry.
        if (allocated(occ_list0)) then
            call enumerate_determinants(.false., .false., sym_in, occ_list0)
        else
            call enumerate_determinants(.false., .false., sym_in)
        end if


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
        ! Don't need to hold determinants, so can just set spawned_size to be 1.
        spawned_size = 1
        allocate(walker_population(1,ndets), stat=ierr)
        call check_allocate('walker_population',ndets,ierr)
        allocate(spawned_walkers1(spawned_size,ndets), stat=ierr)
        call check_allocate('spawned_walkers1',spawned_size*ndets,ierr)
        spawned_walkers => spawned_walkers1
        ! Zero these.
        walker_population = 0
        spawned_walkers = 0

        ! Now we need to set the reference determinant.
        ! We choose the determinant with the lowest Hamiltonian matrix element.
        if (restart) then
            allocate(occ_list0(nel), stat=ierr)
            call check_allocate('occ_list0',nel,ierr)
            allocate(f0(basis_length), stat=ierr)
            call check_allocate('f0',basis_length,ierr)
            call read_restart()
        else
            ref_det = 1
            do i = 2, ndets
                if (hamil(i,i) < hamil(ref_det, ref_det)) then
                    ref_det = i
                end if
            end do

            ! Reference det
            H00 = hamil(ref_det,ref_det)
            if (.not.allocated(f0)) then
                allocate(f0(basis_length), stat=ierr)
                call check_allocate('f0',basis_length,ierr)
            end if
            if (.not.allocated(occ_list0)) then
                allocate(occ_list0(nel), stat=ierr)
                call check_allocate('occ_list0',nel,ierr)
            end if
            call decode_det(f0, occ_list0)
            f0 = dets_list(:,ref_det)
            walker_population(1,ref_det) = nint(D0_population)
        end if

        write (6,'(1X,a29,1X)',advance='no') 'Reference determinant, |D0> ='
        call write_det(dets_list(:,ref_det), new_line=.true.)
        write (6,'(1X,a16,f20.12)') 'E0 = <D0|H|D0> =',H00
        write (6,'(/,1X,a68,/)') 'Note that FCIQMC calculates the correlation energy relative to |D0>.'

    end subroutine init_simple_fciqmc

    subroutine do_simple_fciqmc()

        ! Run the FCIQMC algorithm on the stored Hamiltonian matrix.

        use calc, only: seed
        use determinant_enumeration, only: ndets
        use fciqmc_restart, only: dump_restart, write_restart_file_every_nreports
        use energy_evaluation, only: update_shift
        use parallel, only: parent, iproc
        use utils, only: rng_init_info

        integer :: ireport, icycle, iwalker, ipart
        integer(lint) :: nparticles, nparticles_old
        integer :: nattempts
        real :: t1, t2
        type(dSFMT_t) :: rng

        if (parent) call rng_init_info(seed+iproc)
        call dSFMT_init(seed+iproc, 50000, rng)

        nparticles = sum(abs(walker_population(1,:)))
        nparticles_old = nparticles

        call write_fciqmc_report_header()

        call cpu_time(t1)

        do ireport = 1, nreport

            ! Zero report cycle quantities.
            proj_energy = 0.0_p
            rspawn = 0.0_p
            D0_population = 0.0_p

            do icycle = 1, ncycles

                ! Zero spawning arrays.
                spawned_walkers = 0

                ! Number of spawning attempts that will be made.
                ! convert to integer from integer(lint).  should only be doing
                ! very small calculations with this algothim!
                nattempts = int(nparticles)

                ! Consider all walkers.
                do iwalker = 1, ndets

                    ! It is much easier to evaluate the projected energy at the
                    ! start of the FCIQMC cycle than at the end.
                    call simple_update_proj_energy(iwalker, proj_energy)

                    ! Simulate spawning.
                    do ipart = 1, abs(walker_population(1,iwalker))
                        ! Attempt to spawn from the current particle onto all
                        ! connected determinants.
                        call attempt_spawn(rng, iwalker)
                    end do

                    call simple_death(rng, iwalker)

                end do

                ! Find the spawning rate and add to the running
                ! total.
                rspawn = rspawn + real(sum(abs(spawned_walkers(1,:))))/nattempts

                call simple_annihilation()

            end do

            ! Update the shift
            nparticles = sum(abs(walker_population(1,:)))
            if (vary_shift) then
                call update_shift(nparticles_old, nparticles, ncycles)
            end if
            nparticles_old = nparticles
            if (nparticles > target_particles .and. .not.vary_shift) then
                vary_shift = .true.
            end if

            ! Average these quantities over the report cycle.
            proj_energy = proj_energy/ncycles
            D0_population = D0_population/ncycles
            rspawn = rspawn/ncycles

            call cpu_time(t2)

            ! Output stats
            call write_fciqmc_report(ireport, (/nparticles/), t2-t1, .false.)

            ! Write restart file if required.
            if (mod(ireport,write_restart_file_every_nreports) == 0) &
                call dump_restart(mc_cycles_done+ncycles*ireport, (/nparticles_old/))

            t1 = t2

        end do

        call write_fciqmc_final(ireport)
        write (6,'()')

        if (dump_restart_file) call dump_restart(mc_cycles_done+ncycles*nreport, (/nparticles_old/), vspace=.true.)

    end subroutine do_simple_fciqmc

    subroutine attempt_spawn(rng, iwalker)

        ! Simulate spawning part of FCIQMC algorithm.
        ! We attempt to spawn on all determinants connected to the current
        ! determinant (given by iwalker) with probability tau|K_ij|.  Note this
        ! is different from the optimised FCIQMC algorithm where each walker
        ! only gets one opportunity per FCIQMC cycle to spawn.
        ! In:
        !    iwalker: walker whose particles attempt to clone/die.
        ! In/Out:
        !    rng: random number generator.

        use determinant_enumeration, only: ndets

        integer, intent(in) :: iwalker
        type(dSFMT_t), intent(inout) :: rng

        integer :: j, nspawn
        real(p) :: rate
        real(p) :: r

        ! Simulate spawning by attempting to spawn on all
        ! connected determinants.
        do j = 1, ndets

            ! Can't spawn onto self.
            if (iwalker == j) cycle
            ! Can't spawn onto disconnected dets
            if (hamil(iwalker,j) == 0.0_p) cycle

            ! Attempt spawning.
            ! Spawn with probability tau|K_ij|.
            ! As K_ij = H_ij for off-diagonal elements, we can just use the
            ! stored Hamiltonian matrix directly.
            rate = abs(Tau*hamil(iwalker,j))
            nspawn = int(rate)
            rate = rate - nspawn
            r = get_rand_close_open(rng)
            if (rate > r) nspawn = nspawn + 1

            ! Create particles.
            if (hamil(iwalker,j) > 0.0_p) then
                ! Flip child sign.
                if (walker_population(1,iwalker) < 0) then
                    ! Positive offspring.
                    spawned_walkers(1,j) = spawned_walkers(1,j) + nspawn
                else
                    spawned_walkers(1,j) = spawned_walkers(1,j) - nspawn
                end if
            else
                ! Same sign as parent.
                if (walker_population(1,iwalker) > 0) then
                    ! Positive offspring.
                    spawned_walkers(1,j) = spawned_walkers(1,j) + nspawn
                else
                    spawned_walkers(1,j) = spawned_walkers(1,j) - nspawn
                end if
            end if

        end do

    end subroutine attempt_spawn

    subroutine simple_death(rng, iwalker)

        ! Simulate cloning/death part of FCIQMC algorithm.
        ! In:
        !    iwalker: walker whose particles attempt to clone/die.
        ! In/Out:
        !    rng: random number generator.

        integer, intent(in) :: iwalker
        type(dSFMT_t), intent(inout) :: rng

        integer :: nkill
        real(p) :: rate
        real(dp) :: r

        ! A particle dies with probability, p_d, given by
        !  p_d = tau(K_ii _ S)
        ! where tau is the timestep, S is the shift and K_ii is
        !  K_ii =  < D_i | H | D_i > - E_0
        ! We store the Hamiltonian matrix rather than the K matrix.
        ! It is efficient to allow all particles on a given determinant to
        ! attempt to die in one go (like lemmings) in a stochastic process.
        rate = abs(walker_population(1,iwalker))*tau*(hamil(iwalker,iwalker)-H00-shift)
        ! Number to definitely kill.
        nkill = int(rate)
        rate = rate - nkill

        ! Additional stochasitic death?
        r = get_rand_close_open(rng)
        if (abs(rate) > r) then
            if (rate > 0.0_p) then
                nkill = nkill + 1
            else
                nkill = nkill - 1
            end if
        end if

        ! Don't allow creation of anti-particles in simple_fciqmc.
        if (nkill > abs(walker_population(1,iwalker))) then
            write (6,*) iwalker, walker_population(1,iwalker), &
            abs(walker_population(1,iwalker))*tau*(hamil(iwalker,iwalker)-H00-shift)
            call stop_all('do_simple_fciqmc','Trying to create anti-particles.')
        end if

        ! Update walker populations.
        ! Particle death takes the walker_population closer to 0...
        ! (and similarly if cloning (ie nkill is negative) then the
        ! walker_population should move away from 0...)
        if (walker_population(1,iwalker) > 0) then
            walker_population(1,iwalker) = walker_population(1,iwalker) - nkill
        else
            walker_population(1,iwalker) = walker_population(1,iwalker) + nkill
        end if

    end subroutine simple_death

    subroutine simple_annihilation()

        ! Annihilation: merge main and spawned lists.

        ! This is especially easy as we store the walker populations for all
        ! determinants for both the main and spawned lists so it just amounts to
        ! adding the two arrays together,

        walker_population = walker_population + spawned_walkers

    end subroutine simple_annihilation

    subroutine simple_update_proj_energy(iwalker, inst_proj_energy)

        ! Add the contribution of the current determinant to the projected
        ! energy.
        ! The correlation energy given by the projected energy is:
        !   \sum_{i \neq 0} <D_i|H|D_0> N_i/N_0
        ! where N_i is the population on the i-th determinant, D_i,
        ! and 0 refers to the reference determinant.
        ! During a MC cycle we store
        !   \sum_{i \neq 0} <D_i|H|D_0> N_i
        ! If the current determinant is the reference determinant, then
        ! N_0 is stored as D0_population (defined in fciqmc_data).  This makes
        ! normalisation very efficient.
        ! This procedure is only for the simple fciqmc algorithm, where the
        ! Hamiltonian matrix is explicitly stored.
        ! In:
        !    iwalker: index of current determinant in the main walker list.
        ! In/Out:
        !    inst_proj_energy: running total of the \sum_{i \neq 0} <D_i|H|D_0> N_i.
        !    This is updated if D_i is connected to D_0 (and isn't D_0).

        integer, intent(in) :: iwalker
        real(p), intent(inout) :: inst_proj_energy

        if (iwalker == ref_det) then
            ! Have reference determinant.
            D0_population = D0_population + walker_population(1,iwalker)
        else
            inst_proj_energy = inst_proj_energy + hamil(iwalker,ref_det)*walker_population(1,iwalker)
        end if

    end subroutine simple_update_proj_energy

end module simple_fciqmc
