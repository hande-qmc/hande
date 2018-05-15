program test_mpi

  use, intrinsic :: iso_c_binding
  use mpi_f08

  implicit none

  integer :: origin_count, target_rank, target_count, ierr, rank, msze
  integer :: comm
  integer :: win
  integer :: request
  integer(kind=mpi_address_kind) :: target_disp
  type(c_ptr) :: origin_addr

  call mpi_init(ierr)

  call mpi_comm_size(comm, msze, ierr)
  call mpi_comm_rank(comm, rank, ierr)

  target_rank = 0
  origin_count = 0
  target_disp = 0

  ! This is only for compile time argument checking and does not carry out a
  ! useful or even valid operation
  call mpi_raccumulate(origin_addr, origin_count, mpi_double_precision, target_rank, target_disp, &
                       target_count, mpi_double_precision, mpi_sum, win, request, ierr)

  call mpi_finalize(ierr)

end program
