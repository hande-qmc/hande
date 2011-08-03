module folded_spectrum_utils
!TO DO:
! tau finder
! death step
! spawning step
! calculating P__, P_o, Po_

use const

use proc_pointers

implicit none

! 1)
!      ___
!     / _ \
!    | / \ |
!     \\.//

! 2)
!      __   __ 
!    ./  \./  \.
!
! 3)
!      __
!    ./  \.
!     \__/
!
! 4)
!     _
!    / \
!    \./____.
!
! 5)
!          _
!         / \
!    .____\./


contains



    subroutine fs_spawner(cdet, parent_sign, nspawn, connection)
        ! Attempt to spawn a new particle on a twice connected determinant according to
        ! the fsfciqmc algorithm for a given system
        !
        ! In:
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        !    parent_sign: sign of the population on the parent determinant (i.e.
        !        either a positive or negative integer).
        ! Out:
        !    nspawn: number of particles spawned.  0 indicates the spawning
        !        attempt was unsuccessful.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are spawned.
        use dSFMT_interface, only: genrand_real2
        rng_ptr => genrand_real2 
        
        type(det_info), intent(in) :: cdet
        integer, intent(in) :: parent_sign
        integer, intent(out) :: nspawn(:)
        type(excit), intent(out) :: connection(:)

        real(p),parameter      :: P__=0.2, Po_=(1.0-P__)*0.5, P_o=Po_
        real(p),parameter, dimension(3) :: P_double_elt_type =(/P__, Po_, P_o /)

        real(p)          :: choose_double_elt_type
        
        ! P_gen (k|i) and P_gen (j|k)
        real(p)          :: P_genki, P_genjk 
        real(p)          :: hmatel_ki, hmatel_jk
        real(p)          :: pspawn
        real(p)          :: connection_ki, connection_jk
        real(p)       :: psuccess, pspawn, pgen, hmatel

        ! we need pointers to:
        !      -function that generates an exited determinant for the given system
        !      -function that generates pgen for a given excitation (possibly the same as above
        !      -random number generator function
        !      -

        
        
        ! 1. Choose the type of double element you're going to spawn 
        ! **Do we want a pointer to the random number generating function to determine the type of double element we want to use?**
        choose_double_elt_type = rng_ptr()
        
        ! **We want to choose the largest probability first, since this will reduce the number of if statement calls, 
        ! **however, the values of P__ etc will ideally be reference-determinant dependent, what is the best way to order this sequence?**
        if(choose_double_elt_type <= P__ ) then
            


            
            

        else if (choose_double_elt_type <= P__+P_o ) then
            
            !          _
            !         / \
            !    .____\./
            

            ! 1f. Generate first random excitation and probability of spawning there from cdet 
            call system_excit_ptr(cdet, P_genki, connection_ki)

            ! 2f. Generate the first matrix element.
            call system_matrix_elt_ptr(cdet%f, connection, hmatel_ki)

            ! 1s. Generate the second random excitation 
            !(in this case we stay on the same place)
            P_genjk = 1

            ! 2s. Generate the second matrix element.
            hmatel_jk = !*****stored diagonal matrix element*****
            
            ! 3. Probability of gening...
            pgen = P__ * P_genki * P_genjk
            pspawn = tau*abs(hmatel_ki*hmatel_jk)/pgen
            
            ! 4. Attempt spawning.
            psuccess = rng_ptr()

            ! Need to take into account the possibilty of a spawning attempt
            ! producing multiple offspring...
            ! If pspawn is > 1, then we spawn floor(pspawn) as a minimum and 
            ! then spawn a particle with probability pspawn-floor(pspawn).
            nspawn = int(pspawn)
            pspawn = pspawn - nspawn

            if (pspawn > psuccess) nspawn = nspawn + 1

            if (nspawn > 0) then

                ! 5. If H_ij is positive, then the spawned walker is of opposite
                ! sign to the parent, otherwise the spawned walkers if of the same
                ! sign as the parent.
                if (hmatel_ki*hmatel_jk > 0.0_p) then
                    nspawn = -sign(nspawn, parent_sign)
                else
                    nspawn = sign(nspawn, parent_sign)
                end if

            end if
        else

        endif



    end subroutine fs_spawner





        ! 1. Generate random excitation from cdet and probability of spawning
        ! there.
        call gen_excit_hub_real(cdet, pgen, connection)

        ! 2. find the connecting matrix element.
        call slater_condon1_hub_real_excit(cdet%f, connection, hmatel)

        ! 3. Probability of gening...
        pspawn = tau*abs(hmatel)/pgen

        ! 4. Attempt spawning.
        psuccess = genrand_real2()

        ! Need to take into account the possibilty of a spawning attempt
        ! producing multiple offspring...
        ! If pspawn is > 1, then we spawn floor(pspawn) as a minimum and 
        ! then spawn a particle with probability pspawn-floor(pspawn).
        nspawn = int(pspawn)
        pspawn = pspawn - nspawn

        if (pspawn > psuccess) nspawn = nspawn + 1

        if (nspawn > 0) then

            ! 5. If H_ij is positive, then the spawned walker is of opposite
            ! sign to the parent, otherwise the spawned walkers if of the same
            ! sign as the parent.
            if (hmatel > 0.0_p) then
                nspawn = -sign(nspawn, parent_sign)
            else
                nspawn = sign(nspawn, parent_sign)
            end if

        end if










    subroutine fs_stochastic_death(Kii, population, tot_population, ndeath)

        ! Particles will attempt to die with probability
        !  p_d = tau*M_ii
        ! where tau is the timestep and M_ii is the appropriate diagonal
        ! matrix element.
        ! For folded spectrum FCIQMC, M_ii = K_ii - S where S is the shift and K_ii is
        !  K_ii =  (< D_i | H | D_i > - E_0 - E_offset)^2.
        ! We store < D_i | H | D_i > - E_0 - E_offset, so just need to square it
        ! before passing it to the standard death routine.

        ! In:
        !    Kii: < D_i | H | D_i > - E_0, where D_i is the determinant on
        !         which the particles reside.
        ! In/Out:
        !    population: number of particles on determinant D_i.
        !    tot_population: total number of particles.
        ! Out:
        !    ndeath: running total of number of particles died/cloned.

        use death, only: stochastic_death

        real(p), intent(in) :: Kii
        integer, intent(inout) :: population, tot_population
        integer, intent(out) :: ndeath

        call stochastic_death(Kii**2, population, tot_population, ndeath)

    end subroutine fs_stochastic_death

end module folded_spectrum_utils
