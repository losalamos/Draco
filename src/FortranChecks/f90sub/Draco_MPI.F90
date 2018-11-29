!----------------------------------*-F90-*----------------------------------
!
! \file   wedgehog/f90sub/Draco_MPI.F90
! \author Allan Wollaber
! \date   Mon Jul 30 07:06:24 MDT 2012
! \brief  Helper functions to support scalar vs. distributed MPI tests.
! \note   Copyright (c) 2016-2018 Los Alamos National Security, LLC.
!         All rights reserved.
!
! This is a modified version of jayenne/src/wedgehog/ftest/Wedgehog_MPI.F90
!---------------------------------------------------------------------------

module draco_mpi
  use iso_c_binding, only : c_double, c_intptr_t
  implicit none

  integer, public, save :: f90_rank, f90_num_ranks

  ! Don't pollute the global namespace with the MPI stuff
  private
#ifdef C4_MPI
  include 'mpif.h'
#endif

  ! Make all the subroutines defined below public, though
  public check_mpi_error
  public f90_mpi_init
  public f90_mpi_finalize
  public f90_mpi_barrier

contains

  ! ---------------------------------------------------------------------------
  ! This subroutine checks for and prints the meaning of MPI errors
  ! ---------------------------------------------------------------------------
  subroutine check_mpi_error(error)
    implicit none
    integer, intent(in)             :: error
#ifdef C4_MPI
    integer                         :: error_string_len, ierror
    character(MPI_MAX_ERROR_STRING) :: error_string

    ! Check and report a nonzero error code
    if (error .ne. 0) then
       call mpi_error_string(error, error_string, error_string_len, ierror)
       write(*,"('*** mpi error = ',i18, ' (', a, ')')") error,trim(error_string)
       call MPI_Abort(MPI_COMM_WORLD, 1, ierror)
    end if
#endif
  end subroutine check_mpi_error

  ! ---------------------------------------------------------------------------
  ! A simple MPI initialization function that also sets rank/num_ranks info.
  ! ---------------------------------------------------------------------------
  subroutine f90_mpi_init(ierr)
    implicit none
    integer, intent(out) :: ierr

    ierr = 0
    f90_rank = 0
    f90_num_ranks = 1
#ifdef C4_MPI
    call mpi_init(ierr)
    call check_mpi_error(ierr)

    call mpi_comm_size(MPI_COMM_WORLD, f90_num_ranks, ierr)
    call check_mpi_error(ierr)

    call mpi_comm_rank(MPI_COMM_WORLD, f90_rank, ierr)
    call check_mpi_error(ierr)
#endif
  end subroutine f90_mpi_init

  ! ---------------------------------------------------------------------------
  ! A simple MPI finalize function to wrap the MPI depencencies
  ! ---------------------------------------------------------------------------
  subroutine f90_mpi_finalize(ierr)
    implicit none
    integer, intent(out) :: ierr
#ifdef C4_MPI
    call mpi_finalize(ierr)
    call check_mpi_error(ierr)
#endif
  end subroutine f90_mpi_finalize

  ! ---------------------------------------------------------------------------
  ! Global barrier
  ! ---------------------------------------------------------------------------
  subroutine f90_mpi_barrier(ierr)
    implicit none
    integer, intent(out) :: ierr
#ifdef C4_MPI
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    call check_mpi_error(ierr)
#endif
  end subroutine f90_mpi_barrier

end module draco_mpi

!-----------------------------------------------------------------------------!
! End Draco_MPI.F90
!-----------------------------------------------------------------------------!
