module fciqmc

use fciqmc_data
implicit none

contains

    subroutine init_fciqmc()

        use basis, only: basis_length
        use determinants, only: encode_det, set_spin_polarisation
        use hamiltonian, only: get_hmatel_k
        use system, only: nel
        use hamiltonian, only: slater_condon0_hub_k

        integer :: ierr
        integer :: i, occ_list(nel)

        ! Allocate main walker lists.
        allocate(walker_dets(basis_length,walker_length), stat=ierr)
        allocate(walker_population(walker_length), stat=ierr)
        allocate(walker_energies(walker_length), stat=ierr)

        ! Allocate spawned walker lists.
        allocate(spawned_walker_dets(basis_length,spawned_walker_length), stat=ierr)
        allocate(spawned_walker_population(spawned_walker_length), stat=ierr)

        ! Just for testing...
        ! Ms should become an input option.
        call set_spin_polarisation(0)

        ! Set initial walker population.
        tot_walkers = 1
        walker_population(tot_walkers) = 1
        ! Arbitrary choice initially just for testing!
        forall (i=1:nel) occ_list(i) = i
        walker_dets(:,tot_walkers) = encode_det(occ_list)

        ! Reference determinant.
        walker_energies(tot_walkers) = 0.0_dp

        ! Energy of reference determinant.
        H00 = slater_condon0_hub_k(occ_list)
        
    end subroutine init_fciqmc

    subroutine end_fciqmc()

        ! Deallocate walker lists.
        
        integer :: ierr

        deallocate(walker_dets, stat=ierr)
        deallocate(walker_population, stat=ierr)
        deallocate(walker_energies, stat=ierr)
        deallocate(spawned_walker_dets, stat=ierr)
        deallocate(spawned_walker_population, stat=ierr)

    end subroutine end_fciqmc

    subroutine do_fciqmc()

        use annihilation, only: direct_annihilation
        use basis, only: basis_length
        use death, only: stochastic_death
        use determinants, only: det_info, decode_det_spinocc_spinunocc
        use spawning, only: spawn_hub_k
        use system, only: nel, nalpha, nbeta, nvirt_alpha, nvirt_beta

        integer :: ierr
        integer :: idet, ireport, icycle, iparticles
        type(det_info) :: cdet

        ! Allocate det_info components.
        allocate(cdet%f(basis_length), stat=ierr)
        allocate(cdet%occ_list(nel), stat=ierr)
        allocate(cdet%occ_list_alpha(nalpha), stat=ierr)
        allocate(cdet%occ_list_beta(nbeta), stat=ierr)
        allocate(cdet%unocc_list_alpha(nvirt_alpha), stat=ierr)
        allocate(cdet%unocc_list_beta(nvirt_beta), stat=ierr)

        ! Main fciqmc loop.

        do ireport = 1, nreport

            do icycle = 1, ncycles

                ! Reset the current position in the spawning array to be the
                ! slot preceding the first slot.
                spawning_head = 0
                write (6,*) 'cycle',icycle+(ireport-1)*ncycles

                do idet = 1, tot_walkers ! loop over walkers/dets

                    cdet%f = walker_dets(:,idet)

                    call decode_det_spinocc_spinunocc(cdet%f, cdet)

                    do iparticles = 1, abs(walker_population(idet))
                        
                        ! spawn
                        call spawn_hub_k(cdet, walker_population(idet))

                    end do

                    ! Clone or die.
                    call stochastic_death(idet)

                end do

                call direct_annihilation()

            end do

            ! report

        end do

        write (6,*) 'DONE'

    end subroutine do_fciqmc

end module fciqmc
