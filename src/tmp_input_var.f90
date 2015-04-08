module tmp_input_var

    implicit none

    ! Number of Lanczos eigenpairs to find.
    integer :: nlanczos_eigv = 5

    ! Size of Lanczos basis.
    integer :: lanczos_string_len = 40

    logical :: direct_lanczos = .false.

    logical :: truncate_space = .false.

    ! If true then the determinant list is written to determinant_file.
    logical :: write_determinants = .false.
    character(255) :: determinant_file = 'DETS'

end module tmp_input_var
