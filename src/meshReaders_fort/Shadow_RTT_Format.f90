!----------------------------------*-F90-*----------------------------------
! Shadow_RTT_Format.f90
! B.T. Adams (bta@lanl.gov)
! 27 Sept 99
!---------------------------------------------------------------------------
! @> Shadow_RTT_Format interface file
!---------------------------------------------------------------------------

!===========================================================================
! Shadow_RTT_Format - 
!
! Purpose : Provides shadow interface functions for the C++ RTT_Format class
! for use with Fortran 90 codes. Note that only a destructor is provided
! because the CAR_CU_Interface and a CAR_CU_Build classes (which are shadowed)
! provide full functionality to build an RTT_Format class object, parse an 
! RTT_Format mesh file, and build the mesh.
!
! revision history:
! -----------------
!  0) original
! 
!===========================================================================

      module RTT_Format_Class
          implicit none

          private
          public :: destruct_RTT_Format

!===========================================================================
!
! RTT_Format Class type definition
! 
!===========================================================================

          type, public :: RTT_Format
              integer             :: this
          end type RTT_Format 

!===========================================================================
!
! Interfaces
! 
!===========================================================================

          interface destruct_RTT_Format
              module procedure RTT_Format_destruct
          end interface

          contains

!===========================================================================
!
! Subroutines
! 
!===========================================================================
! Destroy a C++ RTT_Format class object (self).
              subroutine RTT_Format_destruct(self)
                  type(RTT_Format), intent(inout) :: self

                  call destruct_c_rtt_format(self%this)

              end subroutine RTT_Format_destruct

       end module RTT_Format_Class

!---------------------------------------------------------------------------
!                              end of meshReaders_fort/Shadow_RTT_Format.f90
!---------------------------------------------------------------------------
