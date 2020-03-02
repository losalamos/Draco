!----------------------------------*-F90-*----------------------------------//
!!
!* \file   FortranChecks/f90sub/mpi_hw_ftest.f90
!* \author Allan Wollaber, Tom Evans, Kelly Thompson
!* \date   Fri Mar 15 11:25:00 2002
!* \brief  Infinite medium Draco shunt.
!* \note   Copyright (C) 2016-2020 Triad National Security, LLC.
!*         All rights reserved.
!---------------------------------------------------------------------------//

subroutine tst_mpi_hw_f(nf) bind(C, name="tst_mpi_hw")

  use draco_mpi
  use draco_test

  use iso_c_binding, only : c_int, c_double

  implicit none

  integer(c_int), intent(inout) :: nf
  integer :: ierr

  ! ------------------------------------
  ! Initialize MPI
  ! ------------------------------------
  call f90_mpi_init(ierr)
  call check_fail(ierr, f90_rank)

  ! ------------------------------------
  ! Run the problem
  ! ------------------------------------

  if ( f90_rank < f90_num_ranks ) then
     call pass_msg( f90_rank, "MPI rank index ok" )
  else
     call it_fails( f90_rank, "MPI rank > max" )
  endif

  call f90_mpi_barrier(ierr)
  print *, "Hello from rank ", f90_rank, "/", f90_num_ranks

  ! ------------------------------------
  ! Finalize and clean up
  ! ------------------------------------

  call f90_mpi_finalize(ierr)
  call check_fail(ierr, f90_rank)

  ! Print the overall test result
  ! call test_report(f90_rank,nf)

  nf = nf + f90_num_failures

end subroutine tst_mpi_hw_f
