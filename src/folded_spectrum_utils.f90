module folded_spectrum_utils

!TO DO:
! tau finder.
! calculating P__, P_o, Po_.
! optimise diagonal matrix element evaluation.

! Utilities for folded-spectrum method of obtaining excited states in FCIQMC.
! As the Hamiltonian matrix is squared, the spawning/death process must
! stochastically sample both sums in \sum_ij H_ij H_jk.  There are now 4 kinds
! of moves in such a process:

! 1) self spawning
!      ___
!     / _ \
!    | / \ |
!     \\.//

! 2) granddaughter spawning
!      __   __
!    ./  \./  \.
!
!    a special case (which is automatically dealt with) of this is
!      __
!    ./  \.
!     \__/
!
! 3) daughter spawning
!     _
!    / \
!    \./____.
!
! 4) daughter spawning
!          _
!         / \
!    .____\./
!
! FS is easily implemented as a wrapper around FCIQMC for all existing system
! types.
!
! See documentation/theory/folded_spectrum/ for details.

use const

use fciqmc_data, only: fold_line
use determinants, only: det_info_t

implicit none

type(det_info_t), save :: cdet_excit

contains

    subroutine init_folded_spectrum()

        ! Initialise folded spectrum data.

        use fciqmc_data, only: P__, Po_, P_o, X__, Xo_, X_o, tau, H00
        real(p) :: P_renorm

        ! set folded spectrum generation probabilities
        ! renormalise P__, Po_, P_o (just in case)
        P_renorm = P__ + Po_ + P_o
        P__ = P__ / P_renorm
        Po_ = Po_ / P_renorm
        P_o = P_o / P_renorm

        ! calculate chis for split generation
        X__ = sqrt(tau / P__ )
        Xo_ = sqrt(tau / Po_)
        X_o = sqrt(tau / P_o )

        !remove E0 from the fold_line
        fold_line = fold_line - H00

    end subroutine init_folded_spectrum

    subroutine create_cdet_excit(sys, cdet_in, connection, cdet_out)

        ! Generate a complete excited determinant from another determinant and
        ! the excitation information connecting the two determinants.
        !
        ! In:
        !    sys: system being studied.
        !    cdet_in: info on the current determinant that we will excite
        !        from.  The f field must be set.
        !    connection: excitation connecting cdet_in to cdet_out.  Note that
        !        the perm field is not used.
        ! In/Out:
        !    cdet_out info: on the determinant that we will excite to.  Note
        !        that the appropriate fields of cdet_out must be allocated.

        use determinants, only: det_info_t
        use proc_pointers, only: decoder_ptr
        use excitations, only: excit_t, create_excited_det
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        type(det_info_t), intent(in)  :: cdet_in
        type(excit_t), intent(in)     :: connection
        type(det_info_t), intent(inout) :: cdet_out

        ! Create the excited determinant bit string representation
        call create_excited_det(sys%basis, cdet_in%f, connection, cdet_out%f)

        ! Decode the excited determinant bit string representation
        call decoder_ptr(sys, cdet_out%f,cdet_out)

    end subroutine create_cdet_excit

    subroutine fs_spawner(rng, sys, spawn_cutoff, real_factor, cdet, parent_sign, gen_excit_ptr, nspawn, connection)

        ! Attempt to spawn a new particle on a daughter or granddaughter determinant according to
        ! the folded spectrum algorithm for a given system
        !
        ! In/Out:
        !    rng: random number generator.
        ! In:
        !    sys: system being studied.
        !    spawn_cutoff: The size of the minimum spawning event allowed, in
        !        the encoded representation. Events smaller than this will be
        !        stochastically rounded up to this value or down to zero.
        !    real_factor: The factor by which populations are multiplied to
        !        enable non-integer populations.
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        !    parent_sign: sign of the population on the parent determinant (i.e.
        !        either a positive or negative integer).
        !    gen_excit_ptr: procedure pointer to excitation generators.
        !        gen_excit_ptr%full *must* be set to a procedure which generates
        !        a complete excitation.
        ! Out:
        !    nspawn: number of particles spawned, in the encoded representation.
        !        0 indicates the spawning attempt was unsuccessful.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are spawned.

        use excitations, only: excit_t
        use fciqmc_data, only: tau, H00, X__, X_o, Xo_, P__, Po_, P_o
        use excitations, only: create_excited_det, get_excitation
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use system, only: sys_t
        use proc_pointers, only: gen_excit_ptr_t, sc0_ptr
        use spawning, only: nspawn_from_prob, set_child_sign

        implicit none
        type(dSFMT_t), intent(inout) :: rng
        type(sys_t), intent(in) :: sys
        integer(int_p), intent(in) :: spawn_cutoff
        integer(int_p), intent(in) :: real_factor
        type(det_info_t), intent(in) :: cdet
        integer(int_p), intent(in) :: parent_sign
        type(gen_excit_ptr_t), intent(in) :: gen_excit_ptr
        integer(int_p), intent(out) :: nspawn
        type(excit_t), intent(out) :: connection

        real(p)          :: choose_double_elt_type

        real(p)          :: Pgen_ki, Pgen_jk
        real(p)          :: pspawn_ki, pspawn_jk
        integer          :: nspawn_ki, nspawn_jk
        real(p)          :: hmatel_ki, hmatel_jk
        type(excit_t)      :: connection_ki, connection_jk

        ! 0. Choose the type of double element you're going to spawn
        choose_double_elt_type = get_rand_close_open(rng)

        if (choose_double_elt_type <= Po_ ) then
            !     _
            !    / \
            !    \./____.
            !    i,k    j

            ! 1.1 Generate first random excitation and probability of spawning there from cdet
            !    (in this case we stay on the same place)
            Pgen_ki = 1
            hmatel_ki =  sc0_ptr(sys, cdet%f) - H00 - fold_line !***optimise this with stored/calculated values

            ! Calculate P_gen for the first excitation
            pspawn_ki = Xo_ * abs(hmatel_ki) / Pgen_ki

            ! Attempt spawning
            nspawn_ki = nspawn_from_prob(rng, spawn_cutoff, real_factor, pspawn_ki)

            if (nspawn_ki > 0 ) then
            ! Successful spawning on ki

                ! Generate the second random excitation
                call gen_excit_ptr%full(rng, sys, cdet, Pgen_jk, connection_jk, hmatel_jk)

                ! Calculate P_gen for the second excitation
                pspawn_jk = Xo_ * abs(hmatel_jk) / Pgen_jk

                ! Attempt spawning
                nspawn_jk = nspawn_from_prob(rng, spawn_cutoff, real_factor, pspawn_jk)

                if (nspawn_jk > 0) then
                ! Successful spawning on jk

                    ! Combine the spawning numbers
                    nspawn = nspawn_ki * nspawn_jk

                    ! If H_ij is positive, then the spawned walker is of opposite
                    ! sign to the parent, otherwise the spawned walkers if of the same
                    ! sign as the parent.
                    call set_child_sign(hmatel_ki*hmatel_jk, parent_sign, nspawn)

                    ! Set connection to the address of the spawned element
                    connection = connection_jk

                else
                    ! Unsuccessful spawning on jk, set nspawn = 0
                    nspawn = 0
                end if

            else
                ! Unsuccessful spawning  on ki, set nspawn = 0
                nspawn = 0
            end if

        else if (choose_double_elt_type <= Po_ + P_o ) then

            !          _
            !         / \
            !    .____\./
            !    i    k,j

            ! Generate first random excitation and probability of spawning there from cdet
            call gen_excit_ptr%full(rng, sys, cdet, Pgen_ki, connection_ki, hmatel_ki)

            ! Calculate P_gen for the first excitation
            pspawn_ki = X_o * abs(hmatel_ki) / Pgen_ki

            ! Attempt spawning
            nspawn_ki = nspawn_from_prob(rng, spawn_cutoff, real_factor, pspawn_ki)

            if (nspawn_ki > 0 ) then
            ! Successful spawning on ki

                ! Generate the second random excitation
                !    (in this case we stay on the same place)
                ! (i)  generate the first excited determinant
                call create_cdet_excit(sys, cdet, connection_ki, cdet_excit) !could optimise this with create_excited det - we only need %f
                ! (ii) calculate Pgen and hmatel on this site
                Pgen_jk = 1
                hmatel_jk = sc0_ptr(sys, cdet_excit%f) - H00 - fold_line !***optimise this with stored/calculated values

                ! Calculate P_gen for the second excitation
                pspawn_jk = X_o * abs(hmatel_jk) / Pgen_jk

                ! Attempt spawning
                nspawn_jk = nspawn_from_prob(rng, spawn_cutoff, real_factor, pspawn_jk)

                if (nspawn_jk > 0) then
                ! Successful spawning on jk

                    ! Combine the spawning numbers
                    nspawn = nspawn_ki * nspawn_jk

                    ! If H_ij is positive, then the spawned walker is of opposite
                    ! sign to the parent, otherwise the spawned walkers if of the same
                    ! sign as the parent.
                    call set_child_sign(hmatel_ki*hmatel_jk, parent_sign, nspawn)

                    ! Set connection to the address of the spawned element
                    connection = connection_ki

                else
                    ! Unsuccessful spawning on jk, set nspawn = 0
                    nspawn = 0
                end if

            else
                ! Unsuccessful spawning on ki, set nspawn = 0
                nspawn = 0
            end if

        else

            !      __   __
            !    ./  \./  \.
            !    i    k    j
            !
            !    i __
            !    ./  \.
            !     \__/ k
            !    j
            !
            ! The latter is taken care of automatically.

            ! Generate first random excitation and probability of spawning there from cdet
            call gen_excit_ptr%full(rng, sys, cdet, Pgen_ki, connection_ki, hmatel_ki)

            ! Calculate P_gen for the first excitation
            pspawn_ki = X__ * abs(hmatel_ki) / Pgen_ki

            ! Attempt spawning
            nspawn_ki = nspawn_from_prob(rng, spawn_cutoff, real_factor, pspawn_ki)

            if (nspawn_ki > 0 ) then
            ! Successful spawning on ki

                ! Generate the second random excitation
                ! (i)  generate the first excited determinant
                call create_cdet_excit(sys, cdet, connection_ki, cdet_excit)
                ! (ii) excite again
                call gen_excit_ptr%full(rng, sys, cdet_excit, Pgen_jk, connection_jk, hmatel_jk)

                ! Calculate P_gen for the second excitation
                pspawn_jk = X__ * abs(hmatel_jk) / Pgen_jk

                ! Attempt spawning
                nspawn_jk = nspawn_from_prob(rng, spawn_cutoff, real_factor, pspawn_jk)

                if (nspawn_jk > 0) then
                ! Successful spawning on jk

                    ! Combine the spawning numbers
                    nspawn = nspawn_ki * nspawn_jk

                    ! If H_ij is positive, then the spawned walker is of opposite
                    ! sign to the parent, otherwise the spawned walkers if of the same
                    ! sign as the parent.
                    call set_child_sign(hmatel_ki*hmatel_jk, parent_sign, nspawn)

                    ! Calculate the excited determinant (can be up to degree 4)
                    ! (i)   add up the number of excitations
                    connection%nexcit = connection_ki%nexcit + connection_jk%nexcit
                    ! (ii)  combine the annihilations
                    connection%from_orb(:connection_ki%nexcit) = &
                                   connection_ki%from_orb(:connection_ki%nexcit)

                    connection%from_orb(connection_ki%nexcit+1:connection%nexcit) = &
                                   connection_jk%from_orb(:connection_jk%nexcit)

                    ! (iii) combine the creations
                    connection%to_orb(:connection_ki%nexcit) = &
                                   connection_ki%to_orb(:connection_ki%nexcit)

                    connection%to_orb(connection_ki%nexcit+1:connection%nexcit) = &
                                   connection_jk%to_orb(:connection_jk%nexcit)

                    ! (iv) combine the permutations
                    connection%perm = connection_jk%perm .eqv. connection_ki%perm

                else
                    ! Unsuccessful spawning on jk, set nspawn = 0
                    nspawn = 0
                end if

            else
                ! Unsuccessful spawning  on ki, set nspawn = 0
                nspawn = 0
            end if

        endif

    end subroutine fs_spawner

    subroutine fs_stochastic_death(rng, Kii, loc_shift, population, tot_population, ndeath)

        ! Particles will attempt to die with probability
        !  p_d = tau*M_ii
        ! where tau is the timestep and M_ii is the appropriate diagonal
        ! matrix element.
        ! For FSFCIQMC M_ii = (K_ii-fold_line)^2  - S
        ! where S is the shift, fold_line is the point about which we fold
        ! the spectrum and  K_ii is
        !  K_ii =  < D_i | H | D_i > - E_0.

        ! In:
        !    Kii: < D_i | H | D_i > - E_0, where D_i is the determinant on
        !         which the particles reside.
        !    loc_shift: The value of the shift to be used in the death step.
        ! In/Out:
        !    rng: random number generator.
        !    population: number of particles on determinant D_i.
        !    tot_population: total number of particles.
        ! Out:
        !    ndeath: running total of number of particles died/cloned.

        ! Note that population and tot_population refer to a single 'type' of
        ! population, i.e. either a set of Hamiltonian walkers or a set of
        ! Hellmann--Feynman walkers.

        use death, only: stochastic_death
        use dSFMT_interface, only: dSFMT_t

        real(p), intent(in) :: Kii
        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(in) :: loc_shift
        integer(int_p), intent(inout) :: population, ndeath
        real(dp), intent(inout) :: tot_population

        ! Corresponds to:
        !      ___
        !     / _ \
        !    | / \ |
        !     \\.//

        call stochastic_death(rng, (Kii-fold_line)**2, loc_shift, population, tot_population, ndeath)

    end subroutine fs_stochastic_death

end module folded_spectrum_utils
