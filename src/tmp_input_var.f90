module tmp_input_var

    use dmqmc_data, only: rdm_t

    implicit none

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

        ! This stores  information for the various RDMs that the user asks to be
        ! calculated. Each element of this array corresponds to one of these RDMs.
        ! NOTE: This can only be equal to 1 currently.
        type(rdm_t), allocatable :: fci_rdm_info(:)

        ! -- Lanczos only settings --

        ! Number of Lanczos eigenpairs to find.
        integer :: nlanczos_eigv = 5
        ! Size of Lanczos basis.
        integer :: lanczos_string_len = 40
        ! Generate Hamiltonian on the fly (warning: very slow!)
        logical :: direct_lanczos = .false.
    end type fci_in_t

    type(fci_in_t) :: fci_in_global

    logical :: truncate_space = .false.

    ! Use sparse matrix rather than dense matrix?
    logical :: use_sparse_hamil = .false.

end module tmp_input_var
