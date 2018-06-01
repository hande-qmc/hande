module determinant_decoders

! Generation, inspection and manipulation of Slater determinants.

use const
use determinant_data
use parallel, only: parent

implicit none

contains

!--- Decode determinant bit strings ---

    pure subroutine decode_det_occ(sys, f, d, excit_gen_data)

        ! Decode determinant bit string into integer list containing the
        ! occupied orbitals.
        !
        ! In:
        !    sys: system being studied (contains required basis information).
        !    f(tot_string_len): bit string representation of the Slater
        !        determinant.
        !    excit_gen_data (optional): information for excitation generators.
        ! Out:
        !    d: det_info_t variable.  The following components are set:
        !        occ_list: integer list of occupied spin-orbitals in the
        !            Slater determinant.

        use system, only: sys_t
        use excit_gens, only: excit_gen_data_t
        use determinants, only: decode_det

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%tot_string_len)
        type(det_info_t), intent(inout) :: d
        type(excit_gen_data_t), optional, intent(in) :: excit_gen_data

        call decode_det(sys%basis, f, d%occ_list)

    end subroutine decode_det_occ
    
    pure subroutine decode_det_spinocc(sys, f, d, excit_gen_data)

        ! Decode determinant bit string into integer list containing the
        ! occupied orbitals and integer lists containing occupied alpha and
        ! beta orbitals.
        ! In:
        !    sys: system being studied (contains required basis information).
        !    f(tot_string_len): bit string representation of the Slater
        !        determinant.
        !    excit_gen_data (optional): information for excitation generators.
        ! Out:
        !    d: det_info_t variable.  The following components are set:
        !        occ_list: integer list of occupied spin-orbitals in the
        !            Slater determinant.
        !        occ_list_alpha/_beta: integer list of occupied alpha/beta
        !            spin-orbitals in the Slater determinant.

        use system, only: sys_t
        use excit_gens, only: excit_gen_data_t
        use determinants, only: decode_det

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%tot_string_len)
        type(det_info_t), intent(inout) :: d
        type(excit_gen_data_t), optional, intent(in) :: excit_gen_data
        
        integer :: i, ialpha, ibeta

        call decode_det(sys%basis, f, d%occ_list)

        ialpha = 1
        ibeta = 1

        do i = 1, sys%nel
            associate(orb=>d%occ_list(i))
                if (sys%basis%basis_fns(orb)%ms > 0) then
                    d%occ_list_alpha(ialpha) = orb
                    ialpha = ialpha + 1
                else
                    d%occ_list_beta(ibeta) = orb
                    ibeta = ibeta + 1
                end if
            end associate
        end do

    end subroutine decode_det_spinocc
    
    pure subroutine decode_det_occ_unocc(sys, f, d, excit_gen_data)

        ! Decode determinant bit string into integer lists containing the
        ! occupied and unoccupied orbitals.
        !
        ! We return the lists for alpha and beta electrons separately.
        !
        ! In:
        !    sys: system being studied (contains required basis information).
        !    f(tot_string_len): bit string representation of the Slater
        !        determinant.
        !    excit_gen_data (optional): information for excitation generators.
        ! Out:
        !    d: det_info_t variable.  The following components are set:
        !        occ_list: integer list of occupied spin-orbitals in the
        !            Slater determinant.
        !        unocc_list: integer list of unoccupied
        !            spin-orbitals in the Slater determinant.

        use system, only: sys_t
        use excit_gens, only: excit_gen_data_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%tot_string_len)
        type(det_info_t), intent(inout) :: d
        type(excit_gen_data_t), optional, intent(in) :: excit_gen_data
        integer :: i, j, iocc, iunocc, orb, last_basis_ind

        ! A bit too much to do the chunk-based decoding of the occupied list and then fill
        ! in the remaining information.  We only use this in Hubbard model calculations in
        ! k-space, so for now just do a (slow) bit-wise inspection.

        iocc = 0
        iunocc = 0
        orb = 0

        do i = 1, sys%basis%bit_string_len - 1
            ! Manual unrolling allows us to avoid 2 mod statements
            ! and some branching.
            ! [todo] - note above still relevant?
            do j = 0, i0_end
                orb = orb + 1
                if (btest(f(i), j)) then
                    iocc = iocc + 1
                    d%occ_list(iocc) = orb
                else
                    iunocc = iunocc + 1
                    d%unocc_list(iunocc) = orb
                end if
            end do
        end do

        ! Deal with the last element in the determinant bit array separately.
        ! Note that decoding a bit string is surprisingly slow (or, more
        ! importantly, adds up when doing billions of times).
        ! Treating the last element as a special case rather than having an if
        ! statement in the above loop results a speedup of the Hubbard k-space
        ! FCIQMC calculations of 1.5%.
        last_basis_ind = sys%basis%nbasis - i0_length*(sys%basis%bit_string_len-1) - 1
        do j = 0, last_basis_ind
            orb = orb + 1
            if (btest(f(i), j)) then
                iocc = iocc + 1
                d%occ_list(iocc) = orb
            else
                iunocc = iunocc + 1
                d%unocc_list(iunocc) = orb
            end if
        end do

    end subroutine decode_det_occ_unocc

    pure subroutine decode_det_occ_symunocc(sys, f, d, excit_gen_data)

        ! Decode determinant bit string into integer list containing the
        ! occupied orbitals.
        ! In:
        !    sys: system being studied (contains required basis information).
        !    f(tot_string_len): bit string representation of the Slater
        !        determinant.
        !    excit_gen_data (optional): information for excitation generators.
        ! Out:
        !    d: det_info_t variable.  The following components are set:
        !        occ_list: integer list of occupied spin-orbitals in the
        !            Slater determinant.
        !        symunocc(2, sym0_tot:symmax_tot): number of unoccupied orbitals of each
        !            spin/symmetry.  The same indexing scheme is used for
        !            nbasis_sym_spin.

        use system, only: sys_t
        use excit_gens, only: excit_gen_data_t
        use determinants, only: decode_det

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%tot_string_len)
        type(det_info_t), intent(inout) :: d
        type(excit_gen_data_t), optional, intent(in) :: excit_gen_data
        integer :: i, ims, isym

        call decode_det(sys%basis, f, d%occ_list)

        d%symunocc = sys%read_in%pg_sym%nbasis_sym_spin

        do i = 1, sys%nel
            associate(orb=>d%occ_list(i))
                ims = (sys%basis%basis_fns(orb)%ms+3)/2
                isym = sys%basis%basis_fns(orb)%sym
            end associate
            d%symunocc(ims, isym) = d%symunocc(ims, isym) - 1
        end do

    end subroutine decode_det_occ_symunocc

    pure subroutine decode_det_spinocc_symunocc(sys, f, d, excit_gen_data)

        ! Decode determinant bit string into integer list containing the
        ! occupied orbitals.
        ! In:
        !    sys: system being studied (contains required basis information).
        !    f(tot_string_len): bit string representation of the Slater
        !        determinant.
        !    excit_gen_data (optional): information for excitation generators.
        ! Out:
        !    d: det_info_t variable.  The following components are set:
        !        occ_list: integer list of occupied spin-orbitals in the
        !            Slater determinant.
        !        occ_list_alpha/_beta: integer list of occupied alpha/beta
        !            spin-orbitals in the Slater determinant.
        !        symunocc(2, sym0_tot:symmax_tot): number of unoccupied orbitals of each
        !            spin/symmetry.  The same indexing scheme is used for
        !            nbasis_sym_spin.

        use system, only: sys_t
        use excit_gens, only: excit_gen_data_t
        use determinants, only: decode_det

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%tot_string_len)
        type(det_info_t), intent(inout) :: d
        type(excit_gen_data_t), optional, intent(in) :: excit_gen_data
        integer :: i, ims, isym, ialpha, ibeta

        call decode_det(sys%basis, f, d%occ_list)

        d%symunocc = sys%read_in%pg_sym%nbasis_sym_spin
        ialpha = 1
        ibeta = 1

        do i = 1, sys%nel
            associate(orb=>d%occ_list(i))
                ims = (sys%basis%basis_fns(orb)%ms+3)/2
                isym = sys%basis%basis_fns(orb)%sym
                if (ims == 2) then
                    d%occ_list_alpha(ialpha) = orb
                    ialpha = ialpha + 1
                else
                    d%occ_list_beta(ibeta) = orb
                    ibeta = ibeta + 1
                end if
            end associate
            d%symunocc(ims, isym) = d%symunocc(ims, isym) - 1
        end do

    end subroutine decode_det_spinocc_symunocc
    
    pure subroutine decode_det_spinocc_spinunocc(sys, f, d, excit_gen_data)

        ! Decode determinant bit string into integer lists containing the
        ! occupied and unoccupied orbitals.
        !
        ! We return the lists for alpha and beta electrons separately.
        !
        ! In:
        !    sys: system being studied (contains required basis information).
        !    f(tot_string_len): bit string representation of the Slater
        !        determinant.
        !    excit_gen_data (optional): information for excitation generators.
        ! Out:
        !    d: det_info_t variable.  The following components are set:
        !        occ_list: integer list of occupied spin-orbitals in the
        !            Slater determinant.
        !        occ_list_alpha: integer list of occupied alpha
        !            spin-orbitals in the Slater determinant.
        !        occ_list_beta: integer list of occupied beta
        !            spin-orbitals in the Slater determinant.
        !        unocc_list_alpha: integer list of unoccupied alpha
        !            spin-orbitals in the Slater determinant.
        !        unocc_list_beta: integer list of unoccupied beta
        !            spin-orbitals in the Slater determinant.

        use system, only: sys_t
        use excit_gens, only: excit_gen_data_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%tot_string_len)
        type(det_info_t), intent(inout) :: d
        type(excit_gen_data_t), optional, intent(in) :: excit_gen_data
        integer :: i, j, iocc, iocc_a, iocc_b, iunocc_a, iunocc_b, orb, last_basis_ind

        ! A bit too much to do the chunk-based decoding of the occupied list and then fill
        ! in the remaining information.  We only use this in Hubbard model calculations in
        ! k-space, so for now just do a (slow) bit-wise inspection.

        iocc = 0
        iocc_a = 0
        iocc_b = 0
        iunocc_a = 0
        iunocc_b = 0
        orb = 0

        do i = 1, sys%basis%bit_string_len - 1
            ! Manual unrolling allows us to avoid 2 mod statements
            ! and some branching.
            do j = 0, i0_end, 2
                ! Test alpha orbital.
                orb = orb + 1
                if (btest(f(i), j)) then
                    iocc = iocc + 1
                    iocc_a = iocc_a + 1
                    d%occ_list(iocc) = orb
                    d%occ_list_alpha(iocc_a) = orb
                else
                    iunocc_a = iunocc_a + 1
                    d%unocc_list_alpha(iunocc_a) = orb
                end if
                ! Test beta orbital.
                orb = orb + 1
                if (btest(f(i), j+1)) then
                    iocc = iocc + 1
                    iocc_b = iocc_b + 1
                    d%occ_list(iocc) = orb
                    d%occ_list_beta(iocc_b) = orb
                else
                    iunocc_b = iunocc_b + 1
                    d%unocc_list_beta(iunocc_b) = orb
                end if
            end do
        end do

        ! Deal with the last element in the determinant bit array separately.
        ! Note that decoding a bit string is surprisingly slow (or, more
        ! importantly, adds up when doing billions of times).
        ! Treating the last element as a special case rather than having an if
        ! statement in the above loop results a speedup of the Hubbard k-space
        ! FCIQMC calculations of 1.5%.
        last_basis_ind = sys%basis%nbasis - i0_length*(sys%basis%bit_string_len-1) - 1
        do j = 0, last_basis_ind, 2
            ! Test alpha orbital.
            orb = orb + 1
            if (btest(f(i), j)) then
                iocc = iocc + 1
                iocc_a = iocc_a + 1
                d%occ_list(iocc) = orb
                d%occ_list_alpha(iocc_a) = orb
            else
                iunocc_a = iunocc_a + 1
                d%unocc_list_alpha(iunocc_a) = orb
            end if
            ! Test beta orbital.
            orb = orb + 1
            if (btest(f(i), j+1)) then
                iocc = iocc + 1
                iocc_b = iocc_b + 1
                d%occ_list(iocc) = orb
                d%occ_list_beta(iocc_b) = orb
            else
                iunocc_b = iunocc_b + 1
                d%unocc_list_beta(iunocc_b) = orb
            end if
        end do

    end subroutine decode_det_spinocc_spinunocc
    
    pure subroutine decode_det_spinocc_spinsymunocc(sys, f, d, excit_gen_data)

        ! Decode determinant bit string into integer lists containing the
        ! occupied and unoccupied orbitals.
        !
        ! We return the lists for alpha and beta electrons separately.
        !
        ! In:
        !    sys: system being studied (contains required basis information).
        !    f(tot_string_len): bit string representation of the Slater
        !        determinant.
        !    excit_gen_data (optional): information for excitation generators.
        ! Out:
        !    d: det_info_t variable.  The following components are set:
        !        occ_list: integer list of occupied spin-orbitals in the
        !            Slater determinant.
        !        occ_list_alpha: integer list of occupied alpha
        !            spin-orbitals in the Slater determinant.
        !        occ_list_beta: integer list of occupied beta
        !            spin-orbitals in the Slater determinant.
        !        unocc_list_alpha: integer list of unoccupied alpha
        !            spin-orbitals in the Slater determinant.
        !        unocc_list_beta: integer list of unoccupied beta
        !            spin-orbitals in the Slater determinant.
        !        symunocc(2, sym0_tot:symmax_tot): number of unoccupied orbitals of each
        !            spin/symmetry.  The same indexing scheme is used for
        !            nbasis_sym_spin.

        use system, only: sys_t
        use excit_gens, only: excit_gen_data_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%tot_string_len)
        type(det_info_t), intent(inout) :: d
        type(excit_gen_data_t), optional, intent(in) :: excit_gen_data
        integer :: i, j, iocc, iocc_a, iocc_b, iunocc_a, iunocc_b, orb, last_basis_ind, isym, ims

        d%symunocc = sys%read_in%pg_sym%nbasis_sym_spin
        
        ! A bit too much to do the chunk-based decoding of the occupied list and then fill
        ! in the remaining information.

        iocc = 0
        iocc_a = 0
        iocc_b = 0
        iunocc_a = 0
        iunocc_b = 0
        orb = 0

        do i = 1, sys%basis%bit_string_len - 1
            ! Manual unrolling allows us to avoid 2 mod statements
            ! and some branching.
            do j = 0, i0_end, 2
                ! Test alpha orbital.
                orb = orb + 1
                if (btest(f(i), j)) then
                    iocc = iocc + 1
                    iocc_a = iocc_a + 1
                    d%occ_list(iocc) = orb
                    d%occ_list_alpha(iocc_a) = orb
                    ims = 2
                    isym = sys%basis%basis_fns(orb)%sym
                    d%symunocc(ims, isym) = d%symunocc(ims, isym) - 1
                else
                    iunocc_a = iunocc_a + 1
                    d%unocc_list_alpha(iunocc_a) = orb
                end if
                ! Test beta orbital.
                orb = orb + 1
                if (btest(f(i), j+1)) then
                    iocc = iocc + 1
                    iocc_b = iocc_b + 1
                    d%occ_list(iocc) = orb
                    d%occ_list_beta(iocc_b) = orb
                    ims = 1
                    isym = sys%basis%basis_fns(orb)%sym
                    d%symunocc(ims, isym) = d%symunocc(ims, isym) - 1
                else
                    iunocc_b = iunocc_b + 1
                    d%unocc_list_beta(iunocc_b) = orb
                end if
            end do
        end do

        ! Deal with the last element in the determinant bit array separately.
        ! Note that decoding a bit string is surprisingly slow (or, more
        ! importantly, adds up when doing billions of times).
        ! Treating the last element as a special case rather than having an if
        ! statement in the above loop results a speedup of the Hubbard k-space
        ! FCIQMC calculations of 1.5%.
        last_basis_ind = sys%basis%nbasis - i0_length*(sys%basis%bit_string_len-1) - 1
        do j = 0, last_basis_ind, 2
            ! Test alpha orbital.
            orb = orb + 1
            if (btest(f(i), j)) then
                iocc = iocc + 1
                iocc_a = iocc_a + 1
                d%occ_list(iocc) = orb
                d%occ_list_alpha(iocc_a) = orb
                ims = 2
                isym = sys%basis%basis_fns(orb)%sym
                d%symunocc(ims, isym) = d%symunocc(ims, isym) - 1
            else
                iunocc_a = iunocc_a + 1
                d%unocc_list_alpha(iunocc_a) = orb
            end if
            ! Test beta orbital.
            orb = orb + 1
            if (btest(f(i), j+1)) then
                iocc = iocc + 1
                iocc_b = iocc_b + 1
                d%occ_list(iocc) = orb
                d%occ_list_beta(iocc_b) = orb
                ims = 1
                isym = sys%basis%basis_fns(orb)%sym
                d%symunocc(ims, isym) = d%symunocc(ims, isym) - 1
            else
                iunocc_b = iunocc_b + 1
                d%unocc_list_beta(iunocc_b) = orb
            end if
        end do

    end subroutine decode_det_spinocc_spinsymunocc
    
end module determinant_decoders
