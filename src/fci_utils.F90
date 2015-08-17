module fci_utils

! Utility/common procedures and data structures for FCI code.

use dmqmc_data, only: subsys_t

implicit none

! Input options.
type fci_in_t
    ! If true then the non-zero elements of the Hamiltonian matrix are written to hamiltonian_file.
    logical :: write_hamiltonian = .false.
    character(255) :: hamiltonian_file = 'HAMIL'

    ! If true then the determinant list is written to determinant_file.
    logical :: write_determinants = .false.
    character(255) :: determinant_file = 'DETS'

    ! Number of FCI wavefunctions to print out.
    integer :: print_fci_wfn = 0
    ! ...and file to write them to.
    character(255) :: print_fci_wfn_file = 'FCI_WFN'

    ! Number of FCI wavefunctions to compute properties of.
    integer :: analyse_fci_wfn = 0

    ! This stores information for the various RDMs that the user asks to be
    ! calculated. Each element of this array corresponds to one of these RDMs.
    ! NOTE: This can only be equal to 1 currently.
    type(subsys_t), allocatable :: subsys_info(:)

    ! blacs and scalapack split a matrix up into n x n blocks which are then
    ! distributed around the processors in a cyclic fashion.
    ! The block size is critical to performance.  64 seems to be a good value (see
    ! scalapack documentation).
    integer :: block_size = 64

    ! -- Lanczos only settings --

    ! Number of Lanczos eigenpairs to find.
    integer :: nlanczos_eigv = 5
    ! Size of Lanczos basis.
    integer :: lanczos_string_len = 40
    ! Generate Hamiltonian on the fly (warning: very slow!)
    logical :: direct_lanczos = .false.
end type fci_in_t

contains

    subroutine init_fci(sys, fci_in, ref, ndets, dets)

        ! Initialisation of FCI calculations.  (Boilerplate, determinant initialisation.)

        ! In:
        !   sys: system of interest.
        !   fci_in: fci input options.
        ! In/Out:
        !   ref: reference determinant.  If a truncated calculation is being performed and
        !        a reference determinant is not supplied, a simple best guess will be created.
        ! Out:
        !   ndets, dets: number of and bit-string representation of determinants in the
        !        selected Hilbert subspace.

        use const, only: i0

        use determinants, only: spin_orb_list, encode_det, write_det
        use determinant_enumeration, only: enumerate_determinants, print_dets_list
        use qmc_data, only: reference_t
        use reference_determinant, only: set_reference_det
        use symmetry, only: symmetry_orb_list

        use errors, only: stop_all
        use utils, only: get_free_unit, int_fmt
        use parallel, only: parent

        use system, only: sys_t, set_spin_polarisation

        type(sys_t), intent(inout) :: sys
        type(fci_in_t), intent(in) :: fci_in
        type(reference_t), intent(inout) :: ref
        integer, intent(out) :: ndets
        integer(i0), allocatable, intent(out) :: dets(:,:)

        integer, allocatable :: sym_space_size(:)
        integer :: iunit, ref_ms, ref_sym
        logical :: spin_flip
        integer(i0) :: f0(sys%basis%string_len)

        if (parent) then
            write (6,'(1X,"FCI")')
            write (6,'(1X,"---",/)')

            if (fci_in%print_fci_wfn /= 0) then
                ! Overwrite any existing file...
                ! Open a fresh file here so we can just append to it later.
                iunit = get_free_unit()
                open(iunit, file=fci_in%print_fci_wfn_file, status='unknown')
                close(iunit, status='delete')
            end if
        end if

        spin_flip = .false.
        if (allocated(ref%occ_list0)) then
            ! detect if doing spin-flip
            ref_ms = spin_orb_list(sys%basis%basis_fns, ref%occ_list0)
            ref_sym = symmetry_orb_list(sys, ref%occ_list0)
            if (sys%Ms == huge(1)) sys%Ms = ref_ms
            if (sys%symmetry == huge(1)) sys%symmetry = ref_sym
            spin_flip = sys%Ms /= ref_ms
        else if (sys%symmetry == huge(1) .and. sys%nsym == 1) then
            ! Only one option, so don't force it to be set.
            sys%symmetry = sys%sym0
        else if (sys%Ms == huge(1) .or. (sys%symmetry == huge(1))) then
            call stop_all('init_fci', 'Spin and/or symmetry of Hilbert space not defined.')
        end if

        call set_spin_polarisation(sys%basis%nbasis, sys)

        if (parent) call fci_json(sys, fci_in)

        if (.not.allocated(ref%occ_list0) .and. ref%ex_level /= sys%nel) then
            ! Provide a best guess at the reference determinant given symmetry and spin options.
            call set_reference_det(sys, ref%occ_list0, .true., sys%symmetry)
        end if

        ! Construct space
        if (allocated(ref%occ_list0)) then
            call enumerate_determinants(sys, .true., spin_flip, ref%ex_level, sym_space_size, ndets, dets, occ_list0=ref%occ_list0)
            call enumerate_determinants(sys, .false., spin_flip, ref%ex_level, sym_space_size, ndets, dets, &
                                        sys%symmetry, ref%occ_list0)
        else
            call enumerate_determinants(sys, .true., spin_flip, ref%ex_level, sym_space_size, ndets, dets)
            call enumerate_determinants(sys, .false., spin_flip, ref%ex_level, sym_space_size, ndets, dets, sys%symmetry)
        end if
        if (fci_in%write_determinants .and. parent) call print_dets_list(sys, ndets, dets, fci_in%determinant_file)

        ! Info (symmetry, spin, ex_level).
        if (parent) then
            write (6,'(1X,"Symmetry of selected Hilbert subspace:",'//int_fmt(sys%symmetry,1)//',".")') sys%symmetry
            write (6,'(1X,"Spin of selected Hilbert subspace:",'//int_fmt(sys%Ms,1)//',".")') sys%Ms
            if (ref%ex_level /= sys%nel) then
                call encode_det(sys%basis, ref%occ_list0, f0)
                write (6,'(1X,"Reference determinant, |D0> =",1X)',advance='no')
                call write_det(sys%basis, sys%nel, f0, new_line=.true.)
            end if
        end if

    end subroutine init_fci

    subroutine fci_json(sys, fci_in)

        ! Serialise a fci_in_t object in JSON format.

        ! In:
        !   sys: system of interest.
        !   fci_in: fci_in_t object containing fci input values (including any defaults set).

        use json_out
        use dmqmc_data, only: subsys_t_json
        use system, only: sys_t, sys_t_json

        type(sys_t), intent(in) :: sys
        type(fci_in_t), intent(in) :: fci_in
        type(json_out_t) :: js

        call json_object_init(js, tag=.true.)
        call sys_t_json(js, sys)
        call json_object_init(js, 'fci_in')

        call json_write_key(js, 'write_hamiltonian', fci_in%write_hamiltonian)
        call json_write_key(js, 'hamiltonian_file', fci_in%hamiltonian_file)
        call json_write_key(js, 'write_determinants', fci_in%write_determinants)
        call json_write_key(js, 'determinant_file', fci_in%determinant_file)
        call json_write_key(js, 'print_fci_wfn', fci_in%print_fci_wfn)
        call json_write_key(js, 'print_fci_wfn_file', fci_in%print_fci_wfn_file)
        call json_write_key(js, 'analyse_fci_wfn', fci_in%analyse_fci_wfn)
        if (allocated(fci_in%subsys_info)) then
            call subsys_t_json(js, fci_in%subsys_info)
        end if
        call json_write_key(js, 'block_size', fci_in%block_size)
        call json_write_key(js, 'nlanczos_eigv', fci_in%nlanczos_eigv)
        call json_write_key(js, 'lanczos_string_len', fci_in%lanczos_string_len)
        call json_write_key(js, 'direct_lanczos', fci_in%direct_lanczos, .true.)

        call json_object_end(js, terminal=.true.)
        call json_object_end(js, terminal=.true., tag=.true.)
        write (js%io,'()')

    end subroutine fci_json

    subroutine generate_hamil(sys, ndets, dets, hamil, hamil_csr, proc_blacs_info, full_mat)

        ! Generate a symmetry block of the Hamiltonian matrix, H = < D_i | H | D_j >.
        ! The list of determinants, {D_i}, is grouped by symmetry and contains
        ! only determinants of a specified spin.
        ! Only generate the upper diagonal for use with (sca)lapack and Lanczos routines.
        ! In:
        !    sys: system to be studied.
        !    ndets: number of determinants in the Hilbert space.
        !    dets: list of determinants in the Hilbert space (bit-string representation).
        !    proc_blacs_info (optional): BLACS distribution (and related info) of the Hamiltonian matrix.
        !    full_mat (optional): if present and true generate the full matrix rather than
        !        just storing one triangle.
        ! Out:
        !    hamil (optional): Hamiltonian matrix in a square array.
        !    hamil_csr (optional): Hamiltonian matrix in a sparse (compressed sparse row)
        !        format.

        ! Note: either hamil or hamil_csr must be supplied.  If both are supplied then only
        ! hamil is used.

        use checking, only: check_allocate, check_deallocate
        use csr, only: init_csrp, end_csrp, csrp_t
        use utils, only: get_free_unit
        use errors, only: stop_all
        use parallel

        use hamiltonian, only: get_hmatel
        use real_lattice
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: ndets
        integer(i0), intent(in) :: dets(:,:)
        real(p), intent(out), allocatable, optional :: hamil(:,:)
        type(csrp_t), intent(out), optional :: hamil_csr
        type(blacs_info), intent(in), optional :: proc_blacs_info
        logical, intent(in), optional :: full_mat
        integer :: ierr, iunit, n1, n2
        integer :: i, j, ii, jj, ilocal, iglobal, jlocal, jglobal, nnz, imode
        logical :: sparse_mode
        real(p) :: hmatel

        sparse_mode = present(hamil_csr) .and. .not.present(hamil)
        if (.not.present(hamil) .and. .not.present(hamil_csr)) &
            call stop_all('generate_hamil', 'Must supply either hamil or hamil_csr in argument list.')

        if (sparse_mode .and. present(proc_blacs_info)) then
            call stop_all('generate_hamil', &
                'Sparse distributed matrices are not currently implemented.  &
                &If this is disagreeable to you, please contribute patches resolving this situation.')
        end if

        if (present(hamil)) then
            if (present(proc_blacs_info)) then
                allocate(hamil(proc_blacs_info%nrows,proc_blacs_info%ncols), stat=ierr)
            else
                allocate(hamil(ndets, ndets), stat=ierr)
            end if
            call check_allocate('hamil', size(hamil), ierr)
        end if

        ! Form the Hamiltonian matrix < D_i | H | D_j >.
        if (present(proc_blacs_info)) then
            ! blacs divides the matrix up into sub-matrices of size block_size x block_size.
            ! The blocks are distributed in a cyclic fashion amongst the
            ! processors.
            ! Each processor stores a total of nrows of the matrix.
            ! i gives the index of the first row in the current block.
            ! ii gives the index of the current row within the current block.
            ! The local i index is thus i-1+ii.  This is used to refer to the
            ! matrix element as stored (continuously) on the processor.
            ! The global i index is given by the sum of the rows held on
            ! preceeding proecssors for the previous loops over processors (as
            ! done in the loop over i), the rows held on other processors during
            ! the corrent loop over processors and the position within the
            ! current block.
            ! Similarly for the other indices.
            do i = 1, proc_blacs_info%nrows, proc_blacs_info%block_size
                do ii = 1, min(proc_blacs_info%block_size, proc_blacs_info%nrows - i + 1)
                    ilocal = i - 1 + ii
                    iglobal = (i-1)*proc_blacs_info%nproc_rows + proc_blacs_info%procx*proc_blacs_info%block_size + ii
                    do j = 1, proc_blacs_info%ncols, proc_blacs_info%block_size
                        do jj = 1, min(proc_blacs_info%block_size, proc_blacs_info%ncols - j + 1)
                            jlocal = j - 1 + jj
                            jglobal = (j-1)*proc_blacs_info%nproc_cols + proc_blacs_info%procy*proc_blacs_info%block_size + jj
                            hamil(ilocal, jlocal) = get_hmatel(sys, dets(:,iglobal), dets(:,jglobal))
                        end do
                    end do
                end do
            end do
        else if (sparse_mode) then
                ! First, find out how many non-zero elements there are.
                ! We'll be naive and not just test for non-zero by symmetry, but
                ! actually calculate all matrix elements for now.
                ! Then, store the Hamiltonian.
                do imode = 1, 2
                    nnz = 0
                    !$omp parallel
                    do i = 1, ndets
                        ii = i
                        if (present(full_mat)) then
                            if (full_mat) ii = 1
                        end if
                        ! OpenMP chunk size determined completely empirically
                        ! from a single test.  Please feel free to improve...
                        !$omp do private(j, hmatel) schedule(dynamic, 200)
                        do j = ii, ndets
                            hmatel = get_hmatel(sys, dets(:,i), dets(:,j))
                            if (abs(hmatel) > depsilon) then
                                !$omp critical
                                nnz = nnz + 1
                                if (imode == 2) then
                                    hamil_csr%mat(nnz) = hmatel
                                    hamil_csr%col_ind(nnz) = j
                                    if (hamil_csr%row_ptr(i) == 0) hamil_csr%row_ptr(i) = nnz
                                end if
                                !$omp end critical
                            end if
                        end do
                        !$omp end do
                    end do
                    !$omp end parallel
                    if (imode == 1) then
                        write (6,'(1X,a50,i8,/)') 'Number of non-zero elements in Hamiltonian matrix:', nnz
                        call init_csrp(hamil_csr, ndets, nnz, .true.)
                        hamil_csr%row_ptr(1:ndets) = 0
                    else
                        ! Any element not set in row_ptr means that the
                        ! corresponding row has no non-zero elements.
                        ! Set it to be identical to the next row (this avoids
                        ! looping over the zero-row).
                        do i = ndets, 1, -1
                            if (hamil_csr%row_ptr(i) == 0) hamil_csr%row_ptr(i) = hamil_csr%row_ptr(i+1)
                        end do
                    end if
                end do
        else
            !$omp parallel
            do i = 1, ndets
                ii = i
                if (present(full_mat)) then
                    if (full_mat) ii = 1
                end if
                !$omp do private(j) schedule(dynamic, 200)
                do j = ii, ndets
                    hamil(i,j) = get_hmatel(sys,dets(:,i),dets(:,j))
                end do
                !$omp end do
            end do
            !$omp end parallel
        end if

    end subroutine generate_hamil

    subroutine write_hamil(hamiltonian_file, ndets, proc_blacs_info, hamil, hamil_csr)

        ! Write out the Hamiltonian matrix to file.

        ! In:
        !    hamiltonian_file: filename.  Overwritten if exists.
        !    ndets: number of determinants in Hilbert space.
        !    proc_blacs_info (optional, required if nprocs>1): BLACS distribution
        !        (and related info) of the Hamiltonian matrix.
        !    hamil (optional): Hamiltonian matrix in dense matrix format.
        !    hamil_csr (optional): Hamiltonian matrix in sparse (CSR) matrix format.

        use checking, only: check_allocate, check_deallocate
        use const, only: p, depsilon
        use errors, only: stop_all
        use parallel, only: blacs_info, nprocs
        use utils, only: get_free_unit

        use csr, only: csrp_t

        character(*), intent(in) :: hamiltonian_file
        integer, intent(in) :: ndets
        type(blacs_info), intent(in), optional :: proc_blacs_info
        real(p), intent(in), optional :: hamil(:,:)
        type(csrp_t), intent(in), optional :: hamil_csr

        integer :: iunit, i, j, ierr
        real(p), allocatable :: work_print(:)

        iunit = get_free_unit()
        open(iunit, file=hamiltonian_file, status='unknown')
        if (nprocs > 1 .and. present(hamil)) then
            if (.not.present(proc_blacs_info)) call stop_all('write_hamil', 'proc_blacs_info not supplied.')
            ! Note that this uses a different format to the serial case...
            allocate(work_print(proc_blacs_info%block_size**2), stat=ierr)
            call check_allocate('work_print', proc_blacs_info%block_size**2, ierr)
#ifdef PARALLEL
#ifdef SINGLE_PRECISION
            call pslaprnt(ndets, ndets, hamil, 1, 1, proc_blacs_info%desc_m, 0, 0, 'hamil', iunit, work_print)
#else
            call pdlaprnt(ndets, ndets, hamil, 1, 1, proc_blacs_info%desc_m, 0, 0, 'hamil', iunit, work_print)
#endif
#endif
            deallocate(work_print)
            call check_deallocate('work_print', ierr)
        else
            if (present(hamil)) then
                do i=1, ndets
                    write (iunit,*) i,i,hamil(i,i)
                    do j=i+1, ndets
                        if (abs(hamil(i,j)) > depsilon) write (iunit,*) i,j,hamil(i,j)
                    end do
                end do
            else if (present(hamil_csr)) then
                j = 1
                do i = 1, size(hamil_csr%mat)
                    if (abs(hamil_csr%mat(i)) > depsilon) then
                        if (hamil_csr%row_ptr(j+1) <= i) j = j+1
                        write (iunit,*) j, hamil_csr%col_ind(i), hamil_csr%mat(i)
                    end if
                end do
            end if
        end if
        close(iunit, status='keep')

    end subroutine write_hamil

end module fci_utils
