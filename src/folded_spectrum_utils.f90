module folded_spectrum_utils
!TO DO:
! tau finder
! calculating P__, P_o, Po_

use const

use proc_pointers
use fciqmc_data, only: fold_line, fs_offset
use determinants, only: det_info


implicit none

type(det_info), save :: cdet_excit


! 1) self spawning
!      ___
!     / _ \
!    | / \ |
!     \\.//

! 2) granddaughter spawning
!      __   __ 
!    ./  \./  \.
!
! 3) self spawning
!      __
!    ./  \.
!     \__/
!
! 4) daughter spawning
!     _
!    / \
!    \./____.
!
! 5) daughter spawning
!          _
!         / \
!    .____\./


contains

    subroutine init_folded_spectrum()
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

    subroutine alloc_cdet_excit
    use determinants, only: alloc_det_info

        call alloc_det_info(cdet_excit)

    end subroutine alloc_cdet_excit


    subroutine dealloc_cdet_excit
    use determinants, only: dealloc_det_info

        call dealloc_det_info(cdet_excit)

    end subroutine dealloc_cdet_excit

    subroutine create_cdet_excit(cdet_in, connection, cdet_out)
    
        ! Generate a complete excited determinant from another determinant and 
        !the excitation information connecting the two determinants.
        ! In: 
        !    cdet_in: info on the current determinant that we will excite
        !        from.  The f field must be set.
        !    connection: excitation connecting cdet_in to cdet_out.  Note that
        !        the perm field is not used.
        ! Out:
        !    cdet_out info: on the determinant that we will excite to
        use determinants, only: det_info
        use proc_pointers, only: decoder_ptr
        use excitations, only: excit, create_excited_det

        type(det_info), intent(in)  :: cdet_in
        type(excit), intent(in)     :: connection
        type(det_info), intent(inout) :: cdet_out

        ! Create the excited determinant bit string representation
        call create_excited_det(cdet_in%f, connection, cdet_out%f)

        ! Decode the excited determinant bit string representation
        call decoder_ptr(cdet_out%f,cdet_out)

    end subroutine create_cdet_excit










    subroutine fs_spawner(cdet, parent_sign, nspawn, connection)
        ! Attempt to spawn a new particle on a daughter or granddaughter determinant according to
        ! the folded spectrum algorithm for a given system
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
        use excitations, only: excit
        use fciqmc_data, only: tau, H00, X__, X_o, Xo_, P__, Po_, P_o
        use excitations, only: create_excited_det, get_excitation
        use basis, only: basis_length
        use dSFMT_interface, only: genrand_real2

        implicit none
        type(det_info), intent(in) :: cdet
        integer, intent(in) :: parent_sign
        integer, intent(out) :: nspawn
        type(excit), intent(out) :: connection

        real(p)          :: choose_double_elt_type

        real(p)          :: Pgen_ki, Pgen_jk 
        real(p)          :: pspawn_ki, pspawn_jk
        real(p)          :: nspawn_ki, nspawn_jk
        real(p)          :: hmatel_ki, hmatel_jk
        type(excit)      :: connection_ki, connection_jk
        real(p)          :: psuccess

        
        ! 0. Choose the type of double element you're going to spawn 
        choose_double_elt_type = genrand_real2()
        
        if(choose_double_elt_type <= Po_ ) then
            !     _
            !    / \
            !    \./____.
            !    i,k    j

            ! 1.1 Generate first random excitation and probability of spawning there from cdet 
            !    (in this case we stay on the same place)
            Pgen_ki = 1
            hmatel_ki =  sc0_ptr(cdet%f) - H00 - fold_line !***optimise this with stored/calculated values
            ! Calculate P_gen for the first excitation
            pspawn_ki = Xo_ * abs(hmatel_ki) / Pgen_ki

            ! Attempt spawning
            psuccess = genrand_real2()

            ! Multiple offspring
            nspawn_ki = int(pspawn_ki)
            pspawn_ki = pspawn_ki - nspawn_ki

            if (pspawn_ki > psuccess) nspawn_ki = nspawn_ki + 1

            if(nspawn_ki > 0 ) then
            ! Successful spawning on ki

                ! 1.2 Generate the second random excitation 
                call gen_excit_ptr(cdet, Pgen_jk, connection_jk, hmatel_jk)

                ! Calculate P_gen for the second excitation
                pspawn_jk = Xo_ * abs(hmatel_jk) / Pgen_jk
                psuccess = genrand_real2()

                ! Multiple offspring again...
                nspawn_jk = int(pspawn_jk)
                pspawn_jk = pspawn_jk - nspawn_jk


                if (pspawn_jk > psuccess) nspawn_jk = nspawn_jk + 1

                if (nspawn_jk > 0) then
                ! Successful spawning on jk

                    ! Combine the spawning numbers
                    nspawn = nspawn_ki * nspawn_jk

                    ! 4. If H_ij is positive, then the spawned walker is of opposite
                    ! sign to the parent, otherwise the spawned walkers if of the same
                    ! sign as the parent.
                    if (hmatel_ki*hmatel_jk > 0.0_p) then
                        nspawn = -sign(nspawn, parent_sign)
                    else
                        nspawn = sign(nspawn, parent_sign)
                    end if

                    ! 5. Set connection to the address of the spawned element
                    connection = connection_jk

                else
                ! Unsuccessful spawning on jk, set nspawn = 0
                    nspawn =0
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
            call gen_excit_ptr(cdet, Pgen_ki, connection_ki, hmatel_ki)
            

            ! Calculate P_gen for the first excitation
            pspawn_ki = X_o * abs(hmatel_ki) / Pgen_ki

            ! Attempt spawning
            psuccess = genrand_real2()

            ! Multiple offspring
            nspawn_ki = int(pspawn_ki)
            pspawn_ki = pspawn_ki - nspawn_ki

            if (pspawn_ki > psuccess) nspawn_ki = nspawn_ki + 1

            if(nspawn_ki > 0 ) then
            ! Successful spawning on ki

                ! 1.2 Generate the second random excitation 
                !    (in this case we stay on the same place)
                ! (i)  generate the first excited determinant  
                call create_cdet_excit(cdet, connection_ki, cdet_excit) !could optimise this with create_excited det - we only need %f
                ! (ii) calculate Pgen and hmatel on this site       
                Pgen_jk = 1
                hmatel_jk =  sc0_ptr(cdet_excit%f) - H00 - fold_line !***optimise this with stored/calculated values

                ! Calculate P_gen for the second excitation
                pspawn_jk = X_o * abs(hmatel_jk) / Pgen_jk
                psuccess = genrand_real2()

                ! Multiple offspring again...
                nspawn_jk = int(pspawn_jk)
                pspawn_jk = pspawn_jk - nspawn_jk


                if (pspawn_jk > psuccess) nspawn_jk = nspawn_jk + 1

                if (nspawn_jk > 0) then
                ! Successful spawning on jk

                    ! Combine the spawning numbers
                    nspawn = nspawn_ki * nspawn_jk

                    ! 4. If H_ij is positive, then the spawned walker is of opposite
                    ! sign to the parent, otherwise the spawned walkers if of the same
                    ! sign as the parent.
                    if (hmatel_ki*hmatel_jk > 0.0_p) then
                        nspawn = -sign(nspawn, parent_sign)
                    else
                        nspawn = sign(nspawn, parent_sign)
                    end if

                    ! 5. Set connection to the address of the spawned element
                    connection = connection_ki

                else
                ! Unsuccessful spawning on jk, set nspawn = 0
                    nspawn =0
                end if

            else
            ! Unsuccessful spawning  on ki, set nspawn = 0
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

            ! Generate first random excitation and probability of spawning there from cdet 
            call gen_excit_ptr(cdet, Pgen_ki, connection_ki, hmatel_ki)
            

            ! Calculate P_gen for the first excitation
            pspawn_ki = X__ * abs(hmatel_ki) / Pgen_ki

            ! Attempt spawning
            psuccess = genrand_real2()

            ! Multiple offspring
            nspawn_ki = int(pspawn_ki)
            pspawn_ki = pspawn_ki - nspawn_ki

            if (pspawn_ki > psuccess) nspawn_ki = nspawn_ki + 1

            if(nspawn_ki > 0 ) then
            ! Successful spawning on ki

                ! Generate the second random excitation 
                ! (i)  generate the first excited determinant  
                call create_cdet_excit(cdet, connection_ki, cdet_excit)
                ! (ii) excite again
                call gen_excit_ptr(cdet_excit, Pgen_jk, connection_jk, hmatel_jk)

                ! Calculate P_gen for the second excitation
                pspawn_jk = X__ * abs(hmatel_jk) / Pgen_jk
                psuccess = genrand_real2()

                ! Multiple offspring again...
                nspawn_jk = int(pspawn_jk)
                pspawn_jk = pspawn_jk - nspawn_jk


                if (pspawn_jk > psuccess) nspawn_jk = nspawn_jk + 1

                if (nspawn_jk > 0) then
                ! Successful spawning on jk

                    ! Combine the spawning numbers
                    nspawn = nspawn_ki * nspawn_jk

                    ! 4. If H_ij is positive, then the spawned walker is of opposite
                    ! sign to the parent, otherwise the spawned walkers if of the same
                    ! sign as the parent.
                    if (hmatel_ki*hmatel_jk > 0.0_p) then
                        nspawn = -sign(nspawn, parent_sign)
                    else
                        nspawn = sign(nspawn, parent_sign)
                    end if

                ! 5. Calculate the excited determinant (can be up to degree 4)
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
            

                else
                ! Unsuccessful spawning on jk, set nspawn = 0
                    nspawn =0
                end if

            else
            ! Unsuccessful spawning  on ki, set nspawn = 0
                nspawn = 0
            end if



        endif 



    end subroutine fs_spawner










    subroutine fs_stochastic_death(Kii, population, tot_population, ndeath)

        ! Particles will attempt to die with probability
        !  p_d = tau*M_ii
        ! where tau is the timestep and M_ii is the appropriate diagonal
        ! matrix element.
        ! For FSFCIQMC M_ii = (K_ii-fold_line)^2 + fs_offset - S        
        ! where S is the shift, fold_line is the point about which we fold
        ! the spectrum, fs_offset is an addition offset to move the origin 
        ! from zero and  K_ii is
        !  K_ii =  < D_i | H | D_i > - E_0.

        ! In:
        !    Kii: < D_i | H | D_i > - E_0, where D_i is the determinant on
        !         which the particles reside.
        ! In/Out:
        !    population: number of particles on determinant D_i.
        !    tot_population: total number of particles.
        ! Out:
        !    ndeath: running total of number of particles died/cloned.
        
        ! Note that population and tot_population refer to a single 'type' of
        ! population, i.e. either a set of Hamiltonian walkers or a set of
        ! Hellmann--Feynman walkers.
        !      ___
        !     / _ \
        !    | / \ |
        !     \\.//



        use death, only: stochastic_death

        real(p), intent(in) :: Kii
        integer, intent(inout) :: population, tot_population
        integer, intent(out) :: ndeath

        call stochastic_death((Kii-fold_line)**2+fs_offset, population, tot_population, ndeath)

    end subroutine fs_stochastic_death


end module folded_spectrum_utils


