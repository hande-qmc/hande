module tmp_input_var


    use fci_utils, only: fci_in_t

    implicit none

    type(fci_in_t) :: fci_in_global

    logical :: truncate_space = .false.

    ! Use sparse matrix rather than dense matrix?
    logical :: use_sparse_hamil = .false.

    ! Number of report loops use to estimate the energy.
    integer :: nkinetic_cycles

end module tmp_input_var
