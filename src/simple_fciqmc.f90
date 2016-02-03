module simple_fciqmc

! This module contains a very simple, very dumb, self-contained FCIQMC algorithm.
! This is a serial-only algorithm and uses lots of memory---in particular it
! requires that the Hamiltonian matrix and list of determinants are stored.

! Nonetheless, it is useful for debugging and having a simple algorithm which
! definitely works...

use const
use dSFMT_interface
use errors

implicit none

contains

    subroutine init_simple_fciqmc(sys, qmc_in, reference, qs, sparse_hamil, restart, ndets, dets, ref_det, psip_list, spawn, &
                                  hamil)

        ! Initialisation for the simple fciqmc algorithm.
        ! Setup the list of determinants in the space, calculate the relevant
        ! symmetry block of the Hamiltonian matrix, initialise the RNG, allocate
        ! the required memory for the list of walkers and set the initial
        ! walker.

        ! In/Out:
        !    sys: system being studied.  Unaltered on output.
        !    reference: reference determinant. Set on output.
        !    qs: state of QMC calculation.
        ! In:
        !    qmc_in: input options relating to QMC methods.
        !    sparse_hamil: true if using a sparse matrix for the Hamiltonian.
        !    restart: true is restarting from a HDF5 file.
        ! Out:
        !    ndets: number of determinants in the Hilbert space.
        !    dets: list of determinants in the Hilbert space.
        !    ref_det: location of reference in dets.
        !    psip_list: particle_t object for storing psips sampling the Hilbert
        !       space and related information.  Allocated and initialised on
        !       output.
        !    spawn: spawn_t object for holding the spawned psips.  Allocated on
        !       output, with one slot for each determinant in thie Hilbert space.
        !    hamil, hamil_csr: Hamiltonian matrix.  hamil_csr is set if sparse_hamil
        !       is true, otherwise hamil is used.

        use csr, only: csrp_t
        use parallel, only: nprocs
        use checking, only: check_allocate
        use utils, only: int_fmt

        use determinant_enumeration
        use fci_utils, only: generate_hamil, hamil_t
        use qmc_data, only: qmc_in_t, reference_t, particle_t, qmc_state_t
        use spawn_data, only: spawn_t
        use system, only: sys_t, copy_sys_spin_info, set_spin_polarisation

        type(sys_t), intent(inout) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        type(reference_t), intent(inout) :: reference
        type(qmc_state_t), intent(inout) :: qs
        logical, intent(in) :: restart, sparse_hamil
        integer, intent(out) :: ref_det
        type(particle_t), intent(out) :: psip_list
        type(spawn_t), intent(out) :: spawn
        type(hamil_t), intent(out) :: hamil

        integer, allocatable :: sym_space_size(:)
        integer :: ndets
        integer(i0), allocatable :: dets(:,:)

        integer :: ierr
        integer :: i, j
        type(sys_t) :: sys_bak

        if (nprocs > 1) call stop_all('init_simple_fciqmc','Not a parallel algorithm.')

        ! Find and set information about the space.
        call copy_sys_spin_info(sys, sys_bak)
        call set_spin_polarisation(sys%basis%nbasis, sys)
        if (allocated(reference%occ_list0)) then
            call enumerate_determinants(sys, .true., .false., reference%ex_level, sym_space_size, ndets, dets, &
                                        occ_list0=reference%occ_list0)
        else
            call enumerate_determinants(sys, .true., .false., reference%ex_level, sym_space_size, ndets, dets)
        end if

        ! Find all determinants with desired spin and symmetry.
        if (allocated(reference%occ_list0)) then
            call enumerate_determinants(sys, .false., .false., reference%ex_level, sym_space_size, ndets, dets, sys%symmetry, &
                                        reference%occ_list0)
        else
            call enumerate_determinants(sys, .false., .false., reference%ex_level, sym_space_size, ndets, dets, sys%symmetry)
        end if


        ! Set up hamiltonian matrix.
        if (sparse_hamil) then
            call generate_hamil(sys, ndets, dets, hamil, full_mat=.true., use_sparse_hamil=.true.)
        else
            call generate_hamil(sys, ndets, dets, hamil, full_mat=.true.)
        end if

        write (6,'(1X,a13,/,1X,13("-"),/)') 'Simple FCIQMC'
        write (6,'(1X,a53,1X)') 'Using a simple (but correct) serial FCIQMC algorithm.'
        write (6,'(1X,a137)') 'Enumeration of the determinant list and evaluation of &
                              &the Hamiltonian matrix for the given symmetry block and &
                              &spin polarization required.'
        write (6,'(1X,a104,/)') 'This is slow and memory demanding: consider using the &
                                &fciqmc option instead of the simple_fciqmc option.'
        write (6,'(1X,a46,'//int_fmt(sys%symmetry,1)//',1X,a9,'//int_fmt(sys%Ms,1)//',a1,/)') &
            'Considering determinants belonging to symmetry',sys%symmetry,'with spin',sys%Ms,"."

        ! Allocate main and spawned lists to hold population of walkers.
        ! Don't need to hold determinants, so can just set spawned_size to be 1.
        allocate(psip_list%pops(1,ndets), stat=ierr)
        call check_allocate('psip_list%pops',ndets,ierr)
        allocate(spawn%sdata(1,ndets), stat=ierr)
        call check_allocate('spawn%sdata',ndets,ierr)
        ! Zero these.
        psip_list%pops = 0_int_p
        spawn%sdata = 0_int_s

        allocate(qs%shift(1), stat=ierr)
        call check_allocate('qs%shift', 1, ierr)
        qs%shift = qmc_in%initial_shift

        allocate(qs%vary_shift(1), stat=ierr)
        call check_allocate('qs%vary_shift', 1, ierr)
        qs%vary_shift = .false.

        qs%target_particles = qmc_in%target_particles
        qs%tau = qmc_in%tau

        ! Now we need to set the reference determinant.
        ! We choose the determinant with the lowest Hamiltonian matrix element.
        if (restart) then
            allocate(reference%occ_list0(sys%nel), stat=ierr)
            call check_allocate('reference%occ_list0',sys%nel,ierr)
            allocate(reference%f0(sys%basis%string_len), stat=ierr)
            call check_allocate('reference%f0',sys%basis%string_len,ierr)
        else
            if (sparse_hamil) then
                reference%H00 = huge(1.0_p)
                do i = 1, ndets
                    ! mat(k) is M_{ij}, so row_ptr(i) <= k < row_ptr(i+1) and col_ind(k) = j
                    do j = hamil%mat_sparse%row_ptr(i), hamil%mat_sparse%row_ptr(i+1)-1
                        if (hamil%mat_sparse%col_ind(j) == i) then
                            if (hamil%mat_sparse%mat(j) < reference%H00) then
                                reference%H00 = hamil%mat_sparse%mat(j)
                                ref_det = i
                            end if
                        end if
                    end do
                end do
            else
                ref_det = 1
                reference%H00 = hamil%rmat(1,1)
                do i = 2, ndets
                    if (hamil%rmat(i,i) < reference%H00) then
                        ref_det = i
                        reference%H00 = hamil%rmat(ref_det,ref_det)
                    end if
                end do
            end if

            if (.not.allocated(reference%f0)) then
                allocate(reference%f0(sys%basis%string_len), stat=ierr)
                call check_allocate('reference%f0',sys%basis%string_len,ierr)
            end if
            if (.not.allocated(reference%occ_list0)) then
                allocate(reference%occ_list0(sys%nel), stat=ierr)
                call check_allocate('reference%occ_list0',sys%nel,ierr)
            end if
            reference%f0 = dets(:,ref_det)
            call decode_det(sys%basis, reference%f0, reference%occ_list0)
            psip_list%pops(1,ref_det) = nint(qmc_in%D0_population)
        end if

        write (6,'(1X,a29,1X)',advance='no') 'Reference determinant, |D0> ='
        call write_det(sys%basis, sys%nel, dets(:,ref_det), new_line=.true.)
        write (6,'(1X,a16,f20.12)') 'E0 = <D0|H|D0> =',reference%H00
        write (6,'(/,1X,a68,/)') 'Note that FCIQMC calculates the correlation energy relative to |D0>.'

        ! Return sys in an unaltered state.
        call copy_sys_spin_info(sys_bak, sys)

        deallocate(sym_space_size)

    end subroutine init_simple_fciqmc

    subroutine do_simple_fciqmc(sys, qmc_in, restart_in, reference, sparse_hamil, qs, qmc_state_restart)

        ! Run the FCIQMC algorithm on the stored Hamiltonian matrix.

        ! In/Out:
        !    sys: system being studied.  Unaltered on output.
        !    reference: current reference determinant. May be set on input, will be on output.
        ! In:
        !    qmc_in: input options relating to QMC methods.
        !    restart_in: input options for HDF5 restart files.
        !    sparse_hamil: true if using a sparse matrix for the Hamiltonian.
        !    qmc_state_restart (optional): if present, restart from a previous fciqmc calculation
        ! Out:
        !    qs: qmc_state for use if restarting the calculation

        use csr, only: csrp_t
        use json_out

        use energy_evaluation, only: update_shift
        use parallel, only: parent, iproc
        use fciqmc_data, only: write_fciqmc_report_header, write_fciqmc_report
        use qmc_data, only: qmc_in_t, restart_in_t, reference_t, particle_t, qmc_state_t, &
                            qmc_in_t_json, restart_in_t_json, reference_t_json
        use spawn_data, only: spawn_t
        use system, only: sys_t, sys_t_json
        use restart_hdf5, only: dump_restart_hdf5, restart_info_t, init_restart_info_t, dump_restart_file_wrapper
        use check_input, only: check_qmc_opts
        use fci_utils, only: hamil_t

        type(sys_t), intent(inout) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        type(restart_in_t), intent(in) :: restart_in
        type(reference_t), intent(inout) :: reference
        logical, intent(in) :: sparse_hamil
        type(qmc_state_t), intent(out) :: qs
        type(qmc_state_t), intent(in), optional :: qmc_state_restart

        integer :: ireport, icycle, idet, j
        integer(int_p) :: ipart
        real(p) :: nparticles, nparticles_old
        integer :: nattempts
        real :: t1, t2
        type(dSFMT_t) :: rng
        integer :: ref_det, ndets
        integer(i0), allocatable :: dets(:,:)
        real(p) :: H0i, Hii
        type(particle_t) :: psip_list
        type(spawn_t) :: spawn
        logical :: write_restart_shift, restarting
        type(restart_info_t) :: ri, ri_shift
        type(hamil_t) :: hamil
        type(json_out_t) :: js

        ! Check input options
        restarting = present(qmc_state_restart) .or. restart_in%read_restart
        call check_qmc_opts(qmc_in, .false., restarting)

        call init_simple_fciqmc(sys, qmc_in, reference, qs, sparse_hamil, restart_in%read_restart, ndets, dets, ref_det, &
                                psip_list, spawn, hamil)

        if (present(qmc_state_restart)) qs = qmc_state_restart

        if (parent) then
            call json_object_init(js, tag=.true.)
            call sys_t_json(js, sys)
            call qmc_in_t_json(js, qmc_in)
            call restart_in_t_json(js, restart_in)
            call reference_t_json(js, reference, sys)
            call json_write_key(js, 'sparse_hamil', sparse_hamil, .true.)
            call json_object_end(js, terminal=.true., tag=.true.)
            write (js%io,'()')
        end if
        call dSFMT_init(qmc_in%seed+iproc, 50000, rng)

        nparticles = real(sum(abs(psip_list%pops(1,:))),p)
        nparticles_old = nparticles

        call write_fciqmc_report_header(1)

        call cpu_time(t1)
        write_restart_shift = restart_in%write_restart_shift
        call init_restart_info_t(ri, write_id=restart_in%write_id)
        call init_restart_info_t(ri_shift, write_id=restart_in%write_shift_id)

        do ireport = 1, qmc_in%nreport

            ! Zero report cycle quantities.
            qs%estimators%proj_energy = 0.0_p
            qs%spawn_store%rspawn = 0.0_p
            qs%estimators%D0_population = 0.0_p

            do icycle = 1, qmc_in%ncycles

                ! Zero spawning arrays.
                spawn%sdata = 0_int_s

                ! Number of spawning attempts that will be made.
                nattempts = int(nparticles)

                ! Consider all walkers.
                do idet = 1, ndets

                    if (sparse_hamil) then
                        Hii = 0.0_p
                        H0i = 0.0_p
                        do j = hamil%mat_sparse%row_ptr(idet), hamil%mat_sparse%row_ptr(idet+1)-1
                            if (hamil%mat_sparse%col_ind(j) == idet) Hii = hamil%mat_sparse%mat(j)
                            if (hamil%mat_sparse%col_ind(j) == ref_det) H0i = hamil%mat_sparse%mat(j)
                        end do
                    else
                        H0i = hamil%rmat(idet,ref_det)
                        Hii = hamil%rmat(idet,idet)
                    end if

                    ! It is much easier to evaluate the projected energy at the
                    ! start of the FCIQMC cycle than at the end.
                    call simple_update_proj_energy(ref_det == idet, H0i, psip_list%pops(1,idet), qs)

                    ! Attempt to spawn from each particle onto all connected determinants.
                    if (sparse_hamil .and. (hamil%sparse)) then
                        associate(hstart=>hamil%mat_sparse%row_ptr(idet), hend=>hamil%mat_sparse%row_ptr(idet+1)-1)
                            do ipart = 1, abs(psip_list%pops(1,idet))
                                call attempt_spawn(rng, spawn, qs%tau, idet, &
                                                   psip_list%pops(1,idet), hamil%mat_sparse%mat(hstart:hend), &
                                                   hamil%mat_sparse%col_ind(hstart:hend))
                            end do
                        end associate
                    else
                        do ipart = 1, abs(psip_list%pops(1,idet))
                            call attempt_spawn(rng, spawn, qs%tau, idet, psip_list%pops(1,idet), hamil%rmat(:,idet))
                        end do
                    end if

                    call simple_death(rng, qs, Hii, reference%H00, psip_list%pops(1,idet))

                end do

                ! Find the spawning rate and add to the running
                ! total.
                qs%spawn_store%rspawn = qs%spawn_store%rspawn + real(sum(abs(spawn%sdata(1,:))))/nattempts

                call simple_annihilation(spawn, psip_list%pops)

            end do

            ! Update the shift
            nparticles = real(sum(abs(psip_list%pops(1,:))),p)
            if (qs%vary_shift(1)) then
                call update_shift(qmc_in, qs, qs%shift(1), nparticles_old, nparticles, qmc_in%ncycles)
            end if
            nparticles_old = nparticles
            if (nparticles > qs%target_particles .and. .not.qs%vary_shift(1)) then
                qs%vary_shift(1) = .true.
            end if

            ! Average these quantities over the report cycle.
            qs%estimators%proj_energy = qs%estimators%proj_energy/qmc_in%ncycles
            qs%estimators%D0_population = qs%estimators%D0_population/qmc_in%ncycles
            qs%spawn_store%rspawn = qs%spawn_store%rspawn/qmc_in%ncycles

            call cpu_time(t2)

            ! Output stats
            call write_fciqmc_report(qmc_in, qs, ireport, (/nparticles/), t2-t1, .false., .false.)

            ! Write restart file if required.
            call dump_restart_file_wrapper(qs, write_restart_shift, restart_in%write_freq, [nparticles_old], &
                                           ireport, qmc_in%ncycles, sys%basis%nbasis, ri, ri_shift, .false.)

            t1 = t2

        end do

        if (parent) write (6,'()')

        if (restart_in%write_restart) then
            call dump_restart_hdf5(ri, qs, qs%mc_cycles_done+qmc_in%ncycles*qmc_in%nreport, &
                                   (/nparticles_old/), sys%basis%nbasis, .false.)
            if (parent) write (6,'()')
        end if

        deallocate(dets)
        deallocate(spawn%sdata)
        call dSFMT_end(rng)

    end subroutine do_simple_fciqmc

    subroutine attempt_spawn(rng, spawn, tau, idet, pop, hrow, det_indx)

        ! Simulate spawning part of FCIQMC algorithm.
        ! We attempt to spawn on all determinants connected to the current
        ! determinant (given by iwalker) with probability tau|K_ij|.  Note this
        ! is different from the optimised FCIQMC algorithm where each walker
        ! only gets one opportunity per FCIQMC cycle to spawn.

        ! In:
        !    tau: timestep being used.
        !    idet: determinant index from which we attempt to spawn.
        !    pop: population on the parent determinant.
        !    hrow: row of the Hamitlonian matrix i.e. H(:,idet).
        !    det_indx (optional): list of indices for which hrow contains an
        !       entry (i.e. if H(j,idet) is present, then an entry in det_indx
        !       must be set to j.  Must be present if H doesn't contain an entry
        !       for each determinant (ie is in sparse format).
        ! In/Out:
        !    rng: random number generator.
        !    spawn: spawn_t object holding the newly spawned particles.

        use spawn_data, only: spawn_t

        type(dSFMT_t), intent(inout) :: rng
        integer, intent(in) :: idet
        integer(int_p), intent(in) :: pop
        real(p), intent(in) :: tau, hrow(:)
        type(spawn_t), intent(inout) :: spawn
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
            rate = abs(tau*hrow(j))
            nspawn = int(rate, int_s)
            rate = rate - nspawn
            r = get_rand_close_open(rng)
            if (rate > r) nspawn = nspawn + 1_int_s

            ! Create particles.
            if (hrow(j) > 0.0_p) then
                ! Flip child sign.
                if (pop < 0) then
                    ! Positive offspring.
                    spawn%sdata(1,jdet) = spawn%sdata(1,jdet) + nspawn
                else
                    spawn%sdata(1,jdet) = spawn%sdata(1,jdet) - nspawn
                end if
            else
                ! Same sign as parent.
                if (pop > 0) then
                    ! Positive offspring.
                    spawn%sdata(1,jdet) = spawn%sdata(1,jdet) + nspawn
                else
                    spawn%sdata(1,jdet) = spawn%sdata(1,jdet) - nspawn
                end if
            end if

        end do

    end subroutine attempt_spawn

    subroutine simple_death(rng, qs, Hii, H00, pop)

        ! Simulate cloning/death part of FCIQMC algorithm.

        ! In:
        !    Hii: diagonal matrix element, <D_i|H|D_i>
        !    H00: energy of the reference.
        ! In/Out:
        !    rng: random number generator.
        !    pop: population on |D_i>.  On output, the population is updated from applying
        !         the death step.

        use qmc_data, only: qmc_state_t

        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(in) :: Hii, H00
        type(qmc_state_t), intent(in) :: qs
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
        rate = abs(pop)*qs%tau*(Hii-H00-qs%shift(1))
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

    subroutine simple_annihilation(spawn, pop)

        ! Annihilation: merge main and spawned lists.

        ! This is especially easy as we store the populations for all
        ! determinants for both the main and spawned lists so it just amounts to
        ! adding the two arrays together,

        ! In:
        !    spawn: spawn_t object containing the spawned particles.
        ! In/Out:
        !    pop: main population on each site before (input) and after (output)
        !         annihilation.
 
        use spawn_data, only: spawn_t

        type(spawn_t), intent(in)  :: spawn
        integer(int_p), intent(inout) :: pop(:,:)

        pop = pop + int(spawn%sdata, int_p)

    end subroutine simple_annihilation

    subroutine simple_update_proj_energy(ref, H0i, pop, qs)

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
        !    qs: qmc_state_t with running totals of proj_energy, \sum_{i \neq 0} <D_i|H|D_0> N_i,
        !    and the reference population updated.

        use qmc_data, only: qmc_state_t

        logical, intent(in) :: ref
        real(p), intent(in) :: H0i
        integer(int_p), intent(in) :: pop
        type(qmc_state_t), intent(inout) :: qs

        if (ref) then
            ! Have reference determinant.
            qs%estimators%D0_population = qs%estimators%D0_population + pop
        else
            qs%estimators%proj_energy = qs%estimators%proj_energy + H0i*pop
        end if

    end subroutine simple_update_proj_energy

end module simple_fciqmc
