module tmp_input_var

    implicit none

    ! Number of Lanczos eigenpairs to find.
    integer :: nlanczos_eigv = 5

    ! Size of Lanczos basis.
    integer :: lanczos_string_len = 40

    logical :: direct_lanczos = .false.

end module tmp_input_var
