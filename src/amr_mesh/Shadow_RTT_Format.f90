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
! Purpose : Provides shadow interface functions for the C++ Continuous 
! Adaptive Refinement Cartesion Unstructured Mesh RTT_Format class for use
! with Fortran 90 codes. Note that only a destructor is provided because the
! CAR_CU_Interface and a CAR_CU_Build classes (which are shadowed) provide
! full functionality to build an RTT_Format class object, parse an RTT_Format
! mesh file, and build the mesh.
!
! revision history:
! -----------------
!  0) original
! 
!===========================================================================

      module CAR_CU_RTT_Format_Class
          implicit none

          private
          public :: destruct_RTT_Format

!===========================================================================
!
! CAR_CU_RTT_Format Class type definition
! 
!===========================================================================

          type, public :: CAR_CU_RTT_Format
              integer             :: this
          end type CAR_CU_RTT_Format 

!===========================================================================
!
! Interfaces
! 
!===========================================================================

          interface destruct_RTT_Format
              module procedure CAR_CU_RTT_Format_destruct
          end interface

          contains

!===========================================================================
!
! Subroutines
! 
!===========================================================================
! Destroy a C++ CAR_CU_RTT_Format class object (self).
              subroutine CAR_CU_RTT_Format_destruct(self)
                  type(CAR_CU_RTT_Format), intent(inout) :: self

                  call destruct_car_cu_rtt_format(self%this)

              end subroutine CAR_CU_RTT_Format_destruct

       end module CAR_CU_RTT_Format_Class

!---------------------------------------------------------------------------
!                              end of amr_mesh/Shadow_RTT_Format.f90
!---------------------------------------------------------------------------
