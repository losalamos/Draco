!----------------------------------*-F90-*----------------------------------
! Ragged_Right_Array.f90
! B.T. Adams (bta@lanl.gov)
! 29 Oct 99
!---------------------------------------------------------------------------
! @> Ragged_Right_Array derived-type file
!---------------------------------------------------------------------------

!===========================================================================
! Ragged_Right_Array - 
!
! Purpose : Define a F90 derived-type that is equivalent in function to a
!           C++ ragged-right-array.
!
! revision history:
! -----------------
!  0) original
! 
!===========================================================================

      module Ragged_Right_Array
          implicit none

          public 

!===========================================================================
! Class type definitions
!===========================================================================

          type, public :: Integer_Row_Type
              integer, dimension(:), pointer :: column
          end type Integer_Row_Type

          type, public :: Real_Row_Type
              real*8, dimension(:), pointer :: column
          end type Real_Row_Type

          type, public :: Integer_Ragged_Right_Array
              type(Integer_Row_Type), dimension(:), pointer :: row
          end type Integer_Ragged_Right_Array

          type, public :: Real_Ragged_Right_Array
              type(Real_Row_Type), dimension(:), pointer :: row
          end type Real_Ragged_Right_Array

!---------------------------------------------------------------------------
!                              end of Ragged_Right_Array module
!---------------------------------------------------------------------------

      end module Ragged_Right_Array

!---------------------------------------------------------------------------
!                              end of amr_mesh_fort/Ragged_Right_Array.f90
!---------------------------------------------------------------------------
