module symmetry

! Number of symmetries.
! Currently only crystal momentum is implemented.
integer :: nsym

! sym_table(i,j) = k means that k_i + k_j = k_k to within a primitive reciprocal lattice vector.
integer, allocatable :: sym_table(:,:) ! (nsym, nsym)

contains

    subroutine init_symmetry()

        ! Construct the symmetry table.

        use basis, only: nbasis, basis_fns
        use system, only: ndim, system_type, hub_real
        use kpoints, only: is_reciprocal_lattice_vector
        use parallel, only: parent
        use utils, only: int_fmt

        integer :: i, j, ierr
        integer :: ksum(ndim)
        character(10) :: fmt1

        if (system_type == hub_real) then

            nsym = 1

        else

            nsym = nbasis/2
            allocate(sym_table(nsym, nsym), stat=ierr)

            fmt1 = int_fmt(nsym)

            do i = 1, nsym
                do j = i, nsym
                    ksum = basis_fns(i*2)%l + basis_fns(j*2)%l
                    do k = 1, nsym
                        if (is_reciprocal_lattice_vector(ksum - basis_fns(k*2)%l)) then
                            sym_table(i,j) = k
                            sym_table(j,i) = k
                            exit
                        end if
                    end do
                end do
                if (parent) then
                    if (i == 1) then
                        write (6,'(1X,a14,/,1X,14("-"),2/,1X,a82,/,1X,a67,/)') &
                            "Symmetry table", &
                            "The table below gives the result of k_i+k_j to within a reciprocal lattice vector.", &
                            "An index i refers to the wavevector of the i-th alpha spin-orbital."
                    end if
                    do j = 1, nsym
                        write (6,'('//fmt1//')', advance='no') sym_table(j,i)
                    end do
                    write (6,'()')
                    if (i == nsym) write (6,'()')
                end if
            end do

        end if

    end subroutine init_symmetry

    subroutine end_symmetry

        ! Clean up after symmetry.

        integer :: ierr

        if (allocated(sym_table)) deallocate(sym_table, stat=ierr)

    end subroutine end_symmetry

end module symmetry
