module dmqmc_initialise_dm

! This module contains routines to create the starting density matrix, from
! which a DMQMC simulation is performed.

! Currently these routines either sample the infinite-temperature density
! matrix (the identity matrix), or some trial estimate of the density matrix
! at the desired beta value (for the interaction picture DMQMC approach).

! Note that routines for other parts of DMQMC initialisation are kept in
! the dmqmc_procedures module.

use const

implicit none

contains

    subroutine create_initial_density_matrix(rng, sys, qmc_in, dmqmc_in, reference, annihilation_flags, &
                                             target_nparticles_tot, psip_list, spawn)

        ! Create a starting density matrix by sampling the elements of the
        ! (unnormalised) identity matrix. This is a sampling of the
        ! (unnormalised) infinite-temperature density matrix. This is done by
        ! picking determinants/spin configurations with uniform probabilities in
        ! the space being considered.

        ! In/Out:
        !    sys: system being studied. Should be left in an unmodified state on
        !       output
        !    rng: random number generator.
        !    psip_list: particle_t object conaining sample of initial density
        !       matrix on output.
        !    spawn: spawn_t object.  Reset on input and output.  Used to
        !       communicate the generated particles.
        ! In:
        !    qmc_in: input options relating to QMC methods.
        !    dmqmc_in: input options relating to DMQMC.
        !    reference: current reference determinant.
        !    annihilation_flags: calculation specific annihilation flags.
        !    target_nparticles_tot: The total number of psips to attempt to
        !        generate across all processes.

        use annihilation, only: direct_annihilation
        use dSFMT_interface, only:  dSFMT_t, get_rand_close_open
        use errors
        use fciqmc_data, only: real_factor
        use parallel
        use system, only: sys_t, heisenberg, ueg, hub_k, hub_real, read_in, copy_sys_spin_info
        use utils, only: binom_r
        use qmc_common, only: redistribute_particles
        use qmc_data, only: qmc_in_t, reference_t, particle_t, annihilation_flags_t
        use spawn_data, only:spawn_t
        use dmqmc_data, only: dmqmc_in_t

        type(dSFMT_t), intent(inout) :: rng
        type(sys_t), intent(inout) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        type(dmqmc_in_t), intent(in) :: dmqmc_in
        type(reference_t), intent(in) :: reference
        type(annihilation_flags_t), intent(in) :: annihilation_flags
        integer(int_64), intent(in) :: target_nparticles_tot
        type(particle_t), intent(inout) :: psip_list
        type(spawn_t), intent(inout) :: spawn

        real(p) :: nparticles_temp(psip_list%nspaces)
        integer :: nel, ireplica, ierr, ialpha, ms
        integer(int_64) :: npsips_this_proc, npsips
        real(dp) :: total_size, sector_size
        real(dp) :: r, prob
        type(sys_t) :: sys_copy

        npsips_this_proc = target_nparticles_tot/nprocs
        ! If the initial number of psips does not split evenly between all
        ! processors, add the leftover psips to the first processors in order.
        if (target_nparticles_tot-(nprocs*npsips_this_proc) > iproc) &
              npsips_this_proc = npsips_this_proc + 1_int_64

        nparticles_temp = 0.0_p
        total_size = 0

        do ireplica = 1, psip_list%nspaces
            select case(sys%system)
            case(heisenberg)
                if (dmqmc_in%all_spin_sectors) then
                    ! The size (number of configurations) of all symmetry
                    ! sectors combined.
                    total_size = 2.0_dp**(real(sys%lattice%nsites,dp))

                    do nel = 0, sys%lattice%nsites
                        ! The size of this symmetry sector alone.
                        sector_size = binom_r(sys%lattice%nsites, nel)
                        prob = real(npsips_this_proc,dp)*sector_size/total_size
                        npsips = floor(prob, int_64)
                        ! If there are a non-integer number of psips to be
                        ! spawned in this sector then add an extra psip with the
                        ! required probability.
                        prob = prob - npsips
                        r = get_rand_close_open(rng)
                        if (r < prob) npsips = npsips + 1_int_64

                        nparticles_temp(ireplica) = nparticles_temp(ireplica) + real(npsips, p)
                        call random_distribution_heisenberg(rng, sys%basis, nel, npsips, ireplica, spawn)
                    end do
                else
                    ! This process will always create excatly the target number
                    ! of psips.
                    call random_distribution_heisenberg(rng, sys%basis, sys%nel, npsips_this_proc, ireplica, spawn)
                end if
            case(ueg, hub_k, read_in)
                if (dmqmc_in%propagate_to_beta) then
                    ! Initially distribute psips along the diagonal according to
                    ! a guess.
                    if (dmqmc_in%grand_canonical_initialisation) then
                        call init_grand_canonical_ensemble(sys, dmqmc_in, npsips_this_proc, spawn, reference%energy_shift, rng)
                    else
                        call random_distribution_electronic(rng, sys, npsips_this_proc, ireplica, &
                                                                        dmqmc_in%all_sym_sectors, spawn)
                    end if
                    ! Perform metropolis algorithm on initial distribution so
                    ! that we are sampling the trial density matrix.
                    if (dmqmc_in%metropolis_attempts > 0) call initialise_dm_metropolis(sys, rng, qmc_in, dmqmc_in, &
                                                                               npsips_this_proc, ireplica, spawn)
                else
                    if (dmqmc_in%all_spin_sectors) then
                        ! Need to set spin variables appropriately.
                        call copy_sys_spin_info(sys, sys_copy)
                        ! The size (number of configurations) of all spin symmetry
                        ! sectors combined.
                        total_size = 0
                        do ialpha = max(0,sys%nel-sys%basis%nbasis/2), min(sys%nel, sys%basis%nbasis/2)
                            total_size = total_size + &
                                    binom_r(sys%basis%nbasis/2, ialpha)*binom_r(sys%basis%nbasis/2, sys%nel-ialpha)
                        end do

                        do ialpha = max(0,sys%nel-sys%basis%nbasis/2), min(sys%nel, sys%basis%nbasis/2)
                            ! The size of this spin symmetry sector alone.
                            sector_size = binom_r(sys%basis%nbasis/2, ialpha)*binom_r(sys%basis%nbasis/2, sys%nel-ialpha)
                            prob = real(npsips_this_proc,dp)*sector_size/total_size
                            npsips = floor(prob, int_64)
                            ! If there are a non-integer number of psips to be
                            ! spawned in this sector then add an extra psip with the
                            ! required probability.
                            prob = prob - npsips
                            r = get_rand_close_open(rng)
                            if (r < prob) npsips = npsips + 1_int_64

                            nparticles_temp(ireplica) = nparticles_temp(ireplica) + real(npsips, p)
                            ms = 2*ialpha - sys%nel
                            sys%nalpha = (ms + sys%nel) / 2
                            sys%nbeta = sys%nel - sys%nalpha
                            sys%nvirt_alpha = sys%basis%nbasis/2 - sys%nalpha
                            sys%nvirt_beta = sys%basis%nbasis/2 - sys%nbeta
                            call random_distribution_electronic(rng, sys, npsips, ireplica, &
                                                                dmqmc_in%all_sym_sectors, spawn)
                        end do
                        call copy_sys_spin_info(sys_copy, sys)
                    else
                        call random_distribution_electronic(rng, sys, npsips_this_proc, ireplica, &
                                                            dmqmc_in%all_sym_sectors, spawn)
                    end if
                end if
            case(hub_real)
                call random_distribution_electronic(rng, sys, npsips_this_proc, ireplica, dmqmc_in%all_sym_sectors, spawn)
            case default
                call stop_all('create_initial_density_matrix','DMQMC not implemented for this system.')
            end select
        end do

        ! Finally, count the total number of particles across all processes.
        if (dmqmc_in%all_spin_sectors) then
#ifdef PARALLEL
            call mpi_allreduce(nparticles_temp, psip_list%tot_nparticles, psip_list%nspaces, MPI_PREAL, MPI_SUM, &
                                MPI_COMM_WORLD, ierr)
#else
            psip_list%tot_nparticles = nparticles_temp
#endif
        else
            psip_list%tot_nparticles = target_nparticles_tot
        end if

        call direct_annihilation(sys, rng, qmc_in, reference, annihilation_flags, psip_list, spawn)

        if (dmqmc_in%propagate_to_beta) then
            ! Reset the position of the first spawned particle in the spawning array
            spawn%head = spawn%head_start
            ! During the metropolis steps determinants originally in the correct
            ! portions of the spawned walker array are no longer there due to
            ! new determinants being accepted. So we need to reorganise the
            ! determinants appropriately.
            call redistribute_particles(psip_list%states, real_factor, psip_list%pops, psip_list%nstates, &
                                        psip_list%nparticles, spawn)
            call direct_annihilation(sys, rng, qmc_in, reference, annihilation_flags, psip_list, spawn)
        end if

    end subroutine create_initial_density_matrix

    subroutine random_distribution_heisenberg(rng, basis, spins_up, npsips, ireplica, spawn)

        ! For the Heisenberg model only. Distribute the initial number of psips
        ! along the main diagonal. Each diagonal element should be chosen
        ! with the same probability.

        ! Start from a state with all spins down, then choose the above number
        ! of spins to flip up with equal probability.

        ! In/Out:
        !    rng: random number generator.
        !    spawn: spawn_t object to hold spawned particles.
        ! In:
        !    basis: information about the single-particle basis.
        !    spins_up: for the spin configurations generated, this number
        !       specifies how many of the spins shall be up.
        !    npsips: The total number of psips to be created.
        !    ireplica: index of replica (ie which of the possible concurrent
        !       DMQMC populations are we initialising)

        use basis_types, only: basis_t
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use dmqmc_procedures, only: create_diagonal_density_matrix_particle
        use fciqmc_data, only: real_factor
        use spawn_data, only: spawn_t
        use parallel
        use system

        type(dSFMT_t), intent(inout) :: rng
        type(basis_t), intent(in) :: basis
        integer, intent(in) :: spins_up
        integer(int_64), intent(in) :: npsips
        integer, intent(in) :: ireplica
        type(spawn_t), intent(inout) :: spawn
        integer(int_64) :: i
        integer :: rand_basis, bits_set
        integer :: bit_element, bit_position
        integer(i0) :: f(basis%string_len)
        real(dp) :: rand_num

        do i = 1, npsips

            ! Start with all spins down.
            f = 0_i0
            bits_set = 0

            do
                ! If half the spins are now flipped up, we have our basis
                ! function fully created, so exit the loop.
                if (bits_set == spins_up) exit
                ! Choose a random spin to flip.
                rand_num = get_rand_close_open(rng)
                rand_basis = ceiling(rand_num*basis%nbasis)
                ! Find the corresponding positions for this spin.
                bit_position = basis%bit_lookup(1,rand_basis)
                bit_element = basis%bit_lookup(2,rand_basis)
                if (.not. btest(f(bit_element),bit_position)) then
                    ! If not flipped up, flip the spin up.
                    f(bit_element) = ibset(f(bit_element),bit_position)
                    bits_set = bits_set + 1
                end if
            end do

            ! Now call a routine to add the corresponding diagonal element to
            ! the spawned walkers list.
            call create_diagonal_density_matrix_particle(f,basis%string_len, &
                    basis%tensor_label_len, real_factor,ireplica, spawn)

        end do

    end subroutine random_distribution_heisenberg

    subroutine random_distribution_electronic(rng, sys, npsips, ireplica, all_sym_sectors, spawn)

        ! For the electronic Hamiltonians only. Distribute the initial number of psips
        ! along the main diagonal. Each diagonal element should be chosen
        ! with the same probability.

        ! Determinants are generated uniformly in the Hilbert space associated
        ! to selected symmetry sector and spin polarisation.

        ! In:
        !    sys: system being studied.
        !    npsips: The total number of psips to be created.
        !    ireplica: index of replica (ie which of the possible concurrent
        !       DMQMC populations are we initialising)
        !    all_sym_sectors: create determinants in all symmetry sectors?
        ! In/Out:
        !    rng: random number generator
        !    spawn: spawn_t object to hold spawned particles.

        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use symmetry, only: symmetry_orb_list
        use hilbert_space, only: gen_random_det_full_space
        use system, only: sys_t
        use fciqmc_data, only: real_factor
        use spawn_data, only: spawn_t
        use dmqmc_procedures, only: create_diagonal_density_matrix_particle

        type(dSFMT_t), intent(inout) :: rng
        type(sys_t), intent(in) :: sys
        integer(int_64), intent(in) :: npsips
        integer, intent(in) :: ireplica
        logical, intent(in) :: all_sym_sectors
        type(spawn_t), intent(inout) :: spawn

        integer(int_64) :: i
        integer(i0) :: f(sys%basis%string_len)
        integer :: occ_list(sys%nalpha+sys%nbeta)

        do i = 1, npsips
            do
                ! Generate a random determinant uniformly in this specific
                ! symmetry sector and spin polarisation.
                call gen_random_det_full_space(rng, sys, f, occ_list)
                if (all_sym_sectors .or. symmetry_orb_list(sys, occ_list) == sys%symmetry) then
                    call create_diagonal_density_matrix_particle(f, sys%basis%string_len, &
                        sys%basis%tensor_label_len, real_factor, ireplica, spawn)
                    exit
                end if
            end do
        end do

    end subroutine random_distribution_electronic

    subroutine initialise_dm_metropolis(sys, rng, qmc_in, dmqmc_in, npsips, ireplica, spawn)

        ! Attempt to initialise the temperature dependent trial density matrix
        ! using the metropolis algorithm. We either uniformly distribute psips
        ! on all excitation levels or use the grand canonical partition function
        ! as first guess and then use the Metropolis algorithm to distribute
        ! psips according to desired trial density matrix.
        ! The Metropolis algorithm works by generating a new determinant which
        ! is accepted / rejected based on the value of the total energy of the
        ! Slater determinant i.e. E_i = <D_i | H_T | D_i>, where H_T is the
        ! "trial" Hamiltonian.

        ! *** Warning ***: It is up to the user to decide whether enough
        !     metropolis steps have been carried out and that the trial density
        !     matrix is indeed being sampled correctly.

        ! In:
        !    qmc_in: input options relating to QMC methods.
        !    dmqmc_in: input options relating to DMQMC.
        !    npsips: number of psips to distribute in this sector.
        !    ireplica: replica index.
        ! In/Out:
        !    sys: system being studied. Should be left unmodified on output.
        !    rng: random number generator.
        !    spawn: spawn_t object containing the initial distribution of
        !        psips on the diagonal.

        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use system, only: sys_t, copy_sys_spin_info
        use determinants, only: alloc_det_info_t, det_info_t, dealloc_det_info_t, decode_det_spinocc_spinunocc, &
                                encode_det, decode_det, update_sys_spin_info
        use excitations, only: excit_t, create_excited_det
        use fciqmc_data, only: real_factor
        use parallel, only: nprocs, nthreads, parent
        use hilbert_space, only: gen_random_det_truncate_space
        use proc_pointers, only: trial_dm_ptr, gen_excit_ptr, decoder_ptr
        use qmc_data, only: qmc_in_t
        use utils, only: int_fmt
        use spawn_data, only: spawn_t
        use dmqmc_data, only: dmqmc_in_t
        use symmetry, only: symmetry_orb_list

        type(sys_t), intent(inout) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        type(dmqmc_in_t), intent(in) :: dmqmc_in
        integer(int_64), intent(in) :: npsips
        integer, intent(in) :: ireplica
        type(dSFMT_t), intent(inout) :: rng
        type(spawn_t), intent(inout) :: spawn

        integer :: occ_list(sys%nel), naccept
        integer :: idet, iattempt, nsuccess, ms
        integer :: thread_id = 0, proc
        integer(i0) :: f_old(sys%basis%string_len), f_new(sys%basis%string_len)
        real(p), target :: tmp_data(1)
        real(p) :: pgen, hmatel, E_new, E_old, prob
        real(dp) :: r
        type(det_info_t) :: cdet
        type(excit_t) :: connection
        type(sys_t) :: sys_copy

        naccept = 0 ! Number of metropolis moves which are accepted.
        nsuccess = 0 ! Number of successful proposal steps i.e. excluding null excitations.
        idet = 0

        call alloc_det_info_t(sys, cdet)
        ! When averaging over spin the number of virtual orbitals changes so
        ! need to modify them.
        call copy_sys_spin_info(sys, sys_copy)

        ! Visit every psip metropolis_attempts times.
        do iattempt = 1, dmqmc_in%metropolis_attempts
            do proc = 0, nprocs-1
                do idet = spawn%head_start(nthreads-1,proc)+1, spawn%head(thread_id,proc)
                    cdet%f = spawn%sdata(:sys%basis%string_len,idet)
                    E_old = trial_dm_ptr(sys, cdet%f)
                    tmp_data(1) = E_old
                    cdet%data => tmp_data
                    if (dmqmc_in%all_sym_sectors) then
                        call decode_det_spinocc_spinunocc(sys, cdet%f, cdet)
                        if (dmqmc_in%all_spin_sectors) then
                            ! Update spin polarisation properties - these will
                            ! most likely have changed from the previous
                            ! determinant.
                            call update_sys_spin_info(cdet, sys)
                            call dmqmc_spin_flip_metropolis_move(sys, cdet, rng)
                        else
                            call dmqmc_spin_cons_metropolis_move(sys, cdet, rng)
                        end if
                        nsuccess = nsuccess + 1
                        call encode_det(sys%basis, cdet%occ_list, f_new)
                    else
                        call decoder_ptr(sys, cdet%f, cdet)
                        call gen_excit_ptr%full(rng, sys, qmc_in, cdet, pgen, connection, hmatel)
                        ! Check that we didn't generate a null excitation.
                        ! [todo] - Modify accordingly if pgen is ever calculated in for the ueg.
                        if (abs(hmatel) < depsilon) cycle
                        nsuccess = nsuccess + 1
                        call create_excited_det(sys%basis, cdet%f, connection, f_new)
                    end if
                    ! Accept new det with probability p = min[1,exp(-\beta(E_new-E_old))]
                    E_new = trial_dm_ptr(sys, f_new)
                    prob = exp(-1.0_p*dmqmc_in%init_beta*(E_new-E_old))
                    r = get_rand_close_open(rng)
                    if (prob > r) then
                        call decode_det(sys%basis, f_new, occ_list)
                        ! Accept the new determinant by modifying the entry
                        ! in spawned walker list.
                        naccept = naccept + 1
                        spawn%sdata(:sys%basis%string_len,idet) = f_new
                        spawn%sdata(sys%basis%string_len+1:sys%basis%tensor_label_len,idet) = f_new
                    end if
                end do
            end do
        end do

        if (parent) write (6,'(1X,"#",1X, "Average acceptance ratio: ",f8.7,1X," Average number of null excitations: ", f8.7)') &
                           real(naccept)/nsuccess, real(dmqmc_in%metropolis_attempts*npsips-nsuccess)/&
                                                   &(dmqmc_in%metropolis_attempts*npsips)

        call dealloc_det_info_t(cdet)
        call copy_sys_spin_info(sys_copy, sys)

    end subroutine initialise_dm_metropolis

    subroutine dmqmc_spin_flip_metropolis_move(sys, cdet, rng)

        ! Perform Metropolis move where we uniformly select an electron to
        ! excite into an unoccupied orbital while allowing for the total spin of
        ! the determinant to change.

        ! In:
        !    sys: system being studied.
        ! In/Out:
        !    cdet: det_info_t object for the reference determinant.  Must contain
        !        appropriately set occ_list, occ_list_alpha, occ_list_beta,
        !        unocc_list_alpha and unocc_list_beta components. On output
        !        occ_list will be modified.
        !    rng: random number generator.

        use system, only: sys_t
        use determinants, only: det_info_t
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open

        type(sys_t), intent(in) :: sys
        type(det_info_t), intent(inout) :: cdet
        type(dSFMT_t), intent(inout) :: rng

        real(dp) :: species, pflip
        integer :: ielec, iexcit, new_orb

        ! Pick electron at random.
        ielec = int(get_rand_close_open(rng)*sys%nel) + 1
        ! Select unoccupied orbital at random.
        iexcit = int(get_rand_close_open(rng)*sys%nvirt) + 1
        if (iexcit <= sys%nvirt_alpha) then
            new_orb = cdet%unocc_list_alpha(iexcit)
        else
            new_orb = cdet%unocc_list_beta(iexcit-sys%nvirt_alpha)
        end if
        ! Switch the orbitals.
        cdet%occ_list(ielec) = new_orb

    end subroutine dmqmc_spin_flip_metropolis_move

    subroutine dmqmc_spin_cons_metropolis_move(sys, cdet, rng)

        ! Perform Metropolis move where we uniformly select an electron to
        ! excite into an unoccupied orbital while conserving spin.

        ! In:
        !    sys: system being studied.
        ! In/Out:
        !    cdet: det_info_t object for the reference determinant.  Must contain
        !        appropriately set occ_list, occ_list_alpha, occ_list_beta,
        !        unocc_list_alpha and unocc_list_beta components. On output
        !        occ_list will be modified.
        !    rng: random number generator.

        use system, only: sys_t
        use determinants, only: det_info_t
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open

        type(sys_t), intent(in) :: sys
        type(det_info_t), intent(inout) :: cdet
        type(dSFMT_t), intent(inout) :: rng

        real(dp) :: species, pflip
        integer :: orb, ielec, iexcit, new_orb, unocc_list(sys%nvirt)

        ielec = int(get_rand_close_open(rng)*sys%nel) + 1
        orb = cdet%occ_list(ielec)

        ! Conserving ms so only excite inside the correct spin channel as we only
        ! consider single excitations.
        if (mod(orb, 2) == 0) then
            ! Chose a beta spin.
            iexcit = int(get_rand_close_open(rng)*sys%nvirt_beta) + 1
            cdet%occ_list(ielec) = cdet%unocc_list_beta(iexcit)
        else
            ! Chose an alpha spin.
            iexcit = int(get_rand_close_open(rng)*sys%nvirt_alpha) + 1
            cdet%occ_list(ielec) = cdet%unocc_list_alpha(iexcit)
        end if

    end subroutine dmqmc_spin_cons_metropolis_move

    subroutine init_grand_canonical_ensemble(sys, dmqmc_in, npsips, spawn, energy_shift, rng)

        ! Initially distribute psips according to the grand canonical
        ! distribution function.

        ! In:
        !    sys: system being studied.
        !    dmqmc_in: input options for dmqmc.
        !    npsips: number of psips to create on the diagonal.
        !    energy_shift: <D_0|H_new|D_0> - <D_0|H_old|D_0>.
        ! In/Out:
        !    spawn: spawned list.
        !    rng: random number generator.

        use system, only: sys_t
        use spawn_data, only: spawn_t
        use fciqmc_data, only: real_factor
        use symmetry, only: symmetry_orb_list
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use determinants, only: encode_det
        use canonical_kinetic_energy, only: generate_allowed_orbital_list
        use dmqmc_procedures, only: create_diagonal_density_matrix_particle
        use dmqmc_data, only: dmqmc_in_t

        type(sys_t), intent(in) :: sys
        type(dmqmc_in_t), intent(in) :: dmqmc_in
        integer(int_64), intent(in) :: npsips
        real(p), intent(in) :: energy_shift
        type(spawn_t), intent(inout) :: spawn
        type(dSFMT_t), intent(inout) :: rng

        real(dp) :: p_single(sys%basis%nbasis/2)
        real(dp) :: r
        integer :: occ_list(sys%nel)
        integer(i0) :: f(sys%basis%string_len)
        integer :: ireplica, iorb, ipsip
        integer(int_p) :: weight_factor
        logical :: gen

        ireplica = 1

        weight_factor = 1.0

        ! Calculate orbital occupancies.
        ! * Warning *: We assume that we are dealing with a system without
        ! magnetic fields or other funny stuff, so the probabilty of occupying
        ! an alpha spin orbital is equal to that of occupying a beta spin
        ! orbital.
        forall(iorb=1:sys%basis%nbasis:2) p_single(iorb/2+1) = 1.0_p / &
                                          (1+exp(dmqmc_in%init_beta*(sys%basis%basis_fns(iorb)%sp_eigv-sys%ueg%chem_pot)))

        ! In the grand canoical ensemble the probability of occupying a
        ! determinant, |D_i>, is given by \prod_i^N p_i, where the p_i's are the
        ! Fermi factors calculated above in p_single(i). We normally work in the
        ! canonical ensemble, however, so we need to discard any determinant
        ! generated which does not contain nel electrons to obtain the correct
        ! normalisation and hence, the correct distribution. For small systems
        ! the number fluctuations are usually small, so this routine is quite
        ! fast.
        ipsip = 0
        do while (ipsip < npsips)
            occ_list = 0
            ! Select the alpha and beta spin orbitals and discard any
            ! determinant without the correct number of particles.
            if (sys%nalpha > 0) call generate_allowed_orbital_list(sys, rng, p_single, sys%nalpha, 1, occ_list(:sys%nalpha), gen)
            if (.not. gen) cycle
            if (sys%nbeta > 0) call generate_allowed_orbital_list(sys, rng, p_single, sys%nbeta, 0, occ_list(sys%nalpha+1:), gen)
            if (.not. gen) cycle
            ! Create the determinant.
            if (dmqmc_in%all_sym_sectors .or. symmetry_orb_list(sys, occ_list) == sys%symmetry) then
                if (.not. dmqmc_in%free_electron_trial) weight_factor = &
                                                        & calculate_reweighting_factor(sys, occ_list, dmqmc_in%init_beta, energy_shift, rng)
                call encode_det(sys%basis, occ_list, f)
                call create_diagonal_density_matrix_particle(f, sys%basis%string_len, &
                                            sys%basis%tensor_label_len, weight_factor*real_factor, ireplica, spawn)
                ipsip = ipsip + 1
            end if
        end do

        contains

            function calculate_reweighting_factor(sys, occ_list, beta, energy_shift, rng) result(int_weight)

                ! Reweight initial density matrix population so that the desired
                ! diagonal density matrix is sampled.

                ! In:
                !    sys: system being studied.
                !    occ_list: array containing list of occupied orbitals of
                !        generated determinant.
                !    beta: inverse temperature being considered.
                !    energy_shift: <D_0|H_new|D_0> - <D_0|H_old|D_0>.
                ! In/Out:
                !    rng: random number generator.
                ! Returns:
                !    int_weight: integer representation of reweighting
                !        factor.

                use system, only: sys_t
                use proc_pointers, only: energy_diff_ptr
                use dSFMT_interface, only: dSFMT_t, get_rand_close_open

                type(sys_t), intent(in) :: sys
                integer, intent(in) :: occ_list(:)
                real(p), intent(in) :: beta
                real(p), intent(in) :: energy_shift
                type(dSFMT_t), intent(inout) :: rng

                real(p) :: energy_diff, weight, r
                integer(int_p) :: int_weight

                int_weight = 0
                ! Any diagonal density matrix can be sampled from the
                ! non-interacting expression by reweighting, i.e.
                ! p(|D_i>)_new = p(|D_i>)_old * e^{-beta*E_new(|D_i>)}/e^{-beta*E_old(|D_i>)}.
                ! For the case of sampling the "Hartree-Fock" density matrix
                ! this amounts to weighting the populations with the (exponential of the) difference
                ! <D_i|H|D_i>-<D_i|H^0|D_i>, where H^0 is the non-interacting Hamiltonian we use
                ! when doing grand canonical sampling.
                ! We add an (arbitrary) constant energy shift defined above which ensures that p(|D_0>)_new = 1.
                energy_diff = energy_diff_ptr(sys, occ_list)
                weight = exp(-beta*(energy_diff-energy_shift))
                ! Integer part.
                int_weight = int(weight, int_p)
                ! Fractional part.
                weight = weight - int_weight
                r = get_rand_close_open(rng)
                if (weight > r) then
                    int_weight = int_weight + 1
                end if

            end function calculate_reweighting_factor

    end subroutine init_grand_canonical_ensemble

end module dmqmc_initialise_dm
