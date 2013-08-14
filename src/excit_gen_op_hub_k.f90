module excit_gen_op_hub_k

! Module for random excitation generators and related routines for operators
! other than the Hamiltonian for the Hubbard model in momentum space (i.e. Bloch
! orbitals).

! See excit_gen_hub_k for more information about standard excitation generators.
! See hfs_fciqmc for more information about how to sample arbitrary operators
! within FCIQMC.

use const

implicit none

contains

!=== Kinetic energy ===

! The kinetic energy is diagonal in the Bloch basis and so there are no
! determinants connected by the operator to randomly select.

!=== Double occupancy ===

! \hat{D} = 1/L \sum_i n_{i,\uparrow} n_{i,downarrow} (in local orbitals) gives the
! fraction of sites which contain two electrons, where L is the total number of
! sites.  See Becca et al (PRB 61 (2000) R16288).

! In momentum space this becomes (similar to the potential in the Hamiltonian
! operator): 1/L^2 \sum_{k_1,k_2,k_3} c^{\dagger}_{k_1,\uparrow} c^{\dagger}_{k_2,\downarrow} c_{k_3,\downarrow} c_{k_1+k_2-k_3,\uparrow}
! Hence this is trivial to evaluate...it's just like (parts of) the Hamiltonian operator!

! In fact, the connectivity of D is identical to the connectivity of H---the
! only difference is in the value of the matrix elements.  From above, it
! follows that < D | \hat{D} | D' > = < D | \hat{H} | D' > / (U L).

! We can hence abuse the standard excitation generators and the importance
! sampling infrastructure to generate a random excitation and then produce the
! correct matrix element without writing new excitation generators.

    subroutine gen_excit_double_occ_matel_hub_k(cdet, connection, hmatel)

        ! Produce the correct matrix element < D | \hat{D} | D' > for a given
        ! random excitation.

        ! In:
        !    cdet: info on the current determinant, |D>, that we will spawn
        !        from.  Unused; just for interface compatibility.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, |D'>, produced by the excitation
        !        generator.  Unused; just for interface compatibility.
        ! In/Out:
        !    hmatel: on input, the Hamiltonian matrix element, < D | \hat{H} | D'>.
        !        On output, the operator matrix element, < D | \hat{D} | D' >.

        use determinants, only: det_info
        use excitations, only: excit
        use system, only: hubu, nsites

        type(det_info), intent(in) :: cdet
        type(excit), intent(in) :: connection
        real(p), intent(inout) :: hmatel

        hmatel = hmatel / (hubu * nsites)

    end subroutine gen_excit_double_occ_matel_hub_k

end module excit_gen_op_hub_k
