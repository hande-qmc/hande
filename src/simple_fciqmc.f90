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

    subroutine init_simple_fciqmc(sys, ndets, dets, ref_det)

        ! Initialisation for the simple fciqmc algorithm.
        ! Setup the list of determinants in the space, calculate the relevant
        ! symmetry block of the Hamiltonian matrix, initialise the RNG, allocate
        ! the required memory for the list of walkers and set the initial
        ! walker.

        ! In/Out:
        !    sys: system being studied.  Unaltered on output.
        ! Out:
        !    ndets: number of determinants in the Hilbert space.
        !    dets: list of determinants in the Hilbert space.
        !    ref_det: location of reference in dets.

        use parallel, only: nprocs, parent
        use checking, only: check_allocate
        use utils, only: int_fmt

        use determinant_enumeration
        use diagonalisation, only: generate_hamil
        use system, only: sys_t, set_spin_polarisation, copy_sys_spin_info

        type(sys_t), intent(inout) :: sys
        integer, intent(out) :: ref_det

        integer, allocatable :: sym_space_size(:)
        integer :: ndets
        integer(i0), allocatable :: dets(:,:)

        integer :: ierr
        integer :: i, j
        type(sys_t) :: sys_bak

        if (nprocs > 1) call stop_all('init_simple_fciqmc','Not a parallel algorithm.')

        ! Find and set information about the space.
        call copy_sys_spin_info(sys, sys_bak)
        call set_spin_polarisation(sys%basis%nbasis, ms_in, sys)
        if (allocated(occ_list0)) then
            call enumerate_determinants(sys, .true., .false., sym_space_size, ndets, dets, occ_list0=occ_list0)
        else
            call enumerate_determinants(sys, .true., .false., sym_space_size, ndets, dets)
        end if

        ! Find all determinants with desired spin and symmetry.
        if (allocated(occ_list0)) then
            call enumerate_determinants(sys, .false., .false., sym_space_size, ndets, dets, sym_in, occ_list0)
        else
            call enumerate_determinants(sys, .false., .false., sym_space_size, ndets, dets, sym_in)
        end if


        ! Set up hamiltonian matrix.
        call generate_hamil(sys, ndets, dets, use_sparse_hamil, distribute_off)
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
        allocate(walker_population(1,ndets), stat=ierr)
        call check_allocate('walker_population',ndets,ierr)
        allocate(qmc_spawn%sdata(1,ndets), stat=ierr)
        call check_allocate('qmc_spawn%sdata',ndets,ierr)
        ! Zero these.
        walker_population = 0_int_p
        qmc_spawn%sdata = 0_int_s

        allocate(shift(1), stat=ierr)
        call check_allocate('shift', size(shift), ierr)
        shift = initial_shift

        allocate(vary_shift(1), stat=ierr)
        call check_allocate('vary_shift', size(vary_shift), ierr)
        vary_shift = .false.

        ! Now we need to set the reference determinant.
        ! We choose the determinant with the lowest Hamiltonian matrix element.
        if (restart) then
            allocate(occ_list0(sys%nel), stat=ierr)
            call check_allocate('occ_list0',sys%nel,ierr)
            allocate(f0(sys%basis%string_len), stat=ierr)
            call check_allocate('f0',sys%basis%string_len,ierr)
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
                allocate(f0(sys%basis%string_len), stat=ierr)
                call check_allocate('f0',sys%basis%string_len,ierr)
            end if
            if (.not.allocated(occ_list0)) then
                allocate(occ_list0(sys%nel), stat=ierr)
                call check_allocate('occ_list0',sys%nel,ierr)
            end if
            f0 = dets(:,ref_det)
            call decode_det(sys%basis, f0, occ_list0)
            walker_population(1,ref_det) = nint(D0_population)
        end if

        write (6,'(1X,a29,1X)',advance='no') 'Reference determinant, |D0> ='
        call write_det(sys%basis, sys%nel, dets(:,ref_det), new_line=.true.)
        write (6,'(1X,a16,f20.12)') 'E0 = <D0|H|D0> =',H00
        write (6,'(/,1X,a68,/)') 'Note that FCIQMC calculates the correlation energy relative to |D0>.'

        ! Return sys in an unaltered state.
        call copy_sys_spin_info(sys_bak, sys)

        deallocate(sym_space_size)

    end subroutine init_simple_fciqmc

    subroutine do_simple_fciqmc(sys)

        ! Run the FCIQMC algorithm on the stored Hamiltonian matrix.

        ! In/Out:
        !    sys: system being studied.  Unaltered on output.

        use calc, only: seed
        use energy_evaluation, only: update_shift
        use parallel, only: parent, iproc
        use system, only: sys_t
        use utils, only: rng_init_info
        use restart_hdf5, only: dump_restart_hdf5, restart_info_global

        type(sys_t), intent(inout) :: sys
        integer :: ireport, icycle, idet, ipart, j
        real(p) :: nparticles, nparticles_old
        integer :: nattempts
        real :: t1, t2
        type(dSFMT_t) :: rng
        integer :: ref_det, ndets
        integer(i0), allocatable :: dets(:,:)
        real(p) :: H0i, Hii

        call init_simple_fciqmc(sys, ndets, dets, ref_det)

        if (parent) call rng_init_info(seed+iproc)
        call dSFMT_init(seed+iproc, 50000, rng)

        nparticles = real(sum(abs(walker_population(1,:))),p)
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
                qmc_spawn%sdata = 0_int_s

                ! Number of spawning attempts that will be made.
                nattempts = int(nparticles)

                ! Consider all walkers.
                do idet = 1, ndets

                    H0i = hamil(idet,ref_det)
                    Hii = hamil(idet,idet)

                    ! It is much easier to evaluate the projected energy at the
                    ! start of the FCIQMC cycle than at the end.
                    call simple_update_proj_energy(ref_det == idet, H0i, walker_population(1,idet), proj_energy)

                    ! Simulate spawning.
                    do ipart = 1, abs(walker_population(1,idet))
                        ! Attempt to spawn from the current particle onto all
                        ! connected determinants.
                        call attempt_spawn(rng, idet, walker_population(1,idet), hamil(:,idet))
                    end do

                    call simple_death(rng, Hii, walker_population(1,idet))

                end do

                ! Find the spawning rate and add to the running
                ! total.
                rspawn = rspawn + real(sum(abs(qmc_spawn%sdata(1,:))))/nattempts

                call simple_annihilation()

            end do

            ! Update the shift
            nparticles = real(sum(abs(walker_population(1,:))),p)
            if (vary_shift(1)) then
                call update_shift(shift(1), nparticles_old, nparticles, ncycles)
            end if
            nparticles_old = nparticles
            if (nparticles > target_particles .and. .not.vary_shift(1)) then
                vary_shift(1) = .true.
            end if

            ! Average these quantities over the report cycle.
            proj_energy = proj_energy/ncycles
            D0_population = D0_population/ncycles
            rspawn = rspawn/ncycles

            call cpu_time(t2)

            ! Output stats
            call write_fciqmc_report(ireport, (/nparticles/), t2-t1, .false.)

            ! Write restart file if required.
            if (mod(ireport,restart_info_global%write_restart_freq) == 0) &
                call dump_restart_hdf5(restart_info_global, mc_cycles_done+ncycles*ireport, (/nparticles_old/))

            t1 = t2

        end do

        if (parent) write (6,'()')

        if (dump_restart_file) then
            call dump_restart_hdf5(restart_info_global, mc_cycles_done+ncycles*nreport, (/nparticles_old/))
            if (parent) write (6,'()')
        end if

        deallocate(dets)

    end subroutine do_simple_fciqmc

    subroutine attempt_spawn(rng, idet, pop, hrow, det_indx)

        ! Simulate spawning part of FCIQMC algorithm.
        ! We attempt to spawn on all determinants connected to the current
        ! determinant (given by iwalker) with probability tau|K_ij|.  Note this
        ! is different from the optimised FCIQMC algorithm where each walker
        ! only gets one opportunity per FCIQMC cycle to spawn.
        ! In:
        !    iwalker: walker whose particles attempt to clone/die.
        ! In/Out:
        !    rng: random number generator.

        type(dSFMT_t), intent(inout) :: rng
        integer, intent(in) :: idet
        integer(int_p), intent(in) :: pop
        real(p), intent(in) :: hrow(:)
        integer, intent(in), optional :: det_indx(:)

        integer :: j, jdet
        integer(int_s) :: nspawn
        real(p) :: rate
        real(p) :: r

        ! Simulate spawning by attempting to spawn on all
        ! connected determinants.
        do j = 1, ubound(hrow, dim=1)

            if (present(det_indx)) then
                jdet = det_indx(j)
            else
                jdet = j
            end if

            ! Can't spawn onto self.
            if (idet == jdet) cycle
            ! Can't spawn onto disconnected dets
            if (abs(hrow(j)) < depsilon) cycle

            ! Attempt spawning.
            ! Spawn with probability tau|K_ij|.
            ! As K_ij = H_ij for off-diagonal elements, we can just use the
            ! stored Hamiltonian matrix directly.
            rate = abs(Tau*hrow(j))
            nspawn = int(rate, int_s)
            rate = rate - nspawn
            r = get_rand_close_open(rng)
            if (rate > r) nspawn = nspawn + 1_int_s

            ! Create particles.
            if (hrow(j) > 0.0_p) then
                ! Flip child sign.
                if (pop < 0) then
                    ! Positive offspring.
                    qmc_spawn%sdata(1,jdet) = qmc_spawn%sdata(1,jdet) + nspawn
                else
                    qmc_spawn%sdata(1,jdet) = qmc_spawn%sdata(1,jdet) - nspawn
                end if
            else
                ! Same sign as parent.
                if (pop > 0) then
                    ! Positive offspring.
                    qmc_spawn%sdata(1,jdet) = qmc_spawn%sdata(1,jdet) + nspawn
                else
                    qmc_spawn%sdata(1,jdet) = qmc_spawn%sdata(1,jdet) - nspawn
                end if
            end if

        end do

    end subroutine attempt_spawn

    subroutine simple_death(rng, Hii, pop)

        ! Simulate cloning/death part of FCIQMC algorithm.
        ! In:
        !    Hii: diagonal matrix element, <D_i|H|D_i>
        ! In/Out:
        !    rng: random number generator.
        !    pop: population on |D_i>.  On output, the population is updated from applying
        !         the death step.

        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(in) :: Hii
        integer(int_p), intent(inout) :: pop

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
        rate = abs(pop)*tau*(Hii-H00-shift(1))
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
        if (nkill > abs(pop)) then
            write (6,*) pop, abs(pop)*tau*(Hii-H00-shift(1))
            call stop_all('do_simple_fciqmc','Trying to create anti-particles.')
        end if

        ! Update walker populations.
        ! Particle death takes the population closer to 0...
        ! (and similarly if cloning (ie nkill is negative) then the
        ! population should move away from 0...)
        if (pop > 0) then
            pop = pop - nkill
        else
            pop = pop + nkill
        end if

    end subroutine simple_death

    subroutine simple_annihilation()

        ! Annihilation: merge main and spawned lists.

        ! This is especially easy as we store the walker populations for all
        ! determinants for both the main and spawned lists so it just amounts to
        ! adding the two arrays together,

        walker_population = walker_population + int(qmc_spawn%sdata, int_p)

    end subroutine simple_annihilation

    subroutine simple_update_proj_energy(ref, H0i, pop, proj_energy)

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
        !    ref: true if |D_i> is the reference, |D_0>.
        !    pop: population on |D_i>.
        ! In/Out:
        !    proj_energy: running total of the \sum_{i \neq 0} <D_i|H|D_0> N_i.
        !    This is updated if |D_i> is connected to |D_0> (and isn't |D_0>).

        logical, intent(in) :: ref
        real(p), intent(in) :: H0i
        integer(int_p), intent(in) :: pop
        real(p), intent(inout) :: proj_energy

        if (ref) then
            ! Have reference determinant.
            D0_population = D0_population + pop
        else
            proj_energy = proj_energy + H0i*pop
        end if

    end subroutine simple_update_proj_energy

end module simple_fciqmc
