!----------------------------------*-F90-*----------------------------------
! Shadow_CAR_CU_Interface.f90
! B.T. Adams (bta@lanl.gov)
! 27 Sept 99
!---------------------------------------------------------------------------
! @> Shadow_CAR_CU_Interface interface file
!---------------------------------------------------------------------------

!===========================================================================
! Shadow_CAR_CU_Interface - 
!
! Purpose : Provides shadow interface functions for the C++ Continuous 
! Adaptive Refinement Cartesion Unstructured Mesh Interface Class
!
! revision history:
! -----------------
!  0) original
! 
!===========================================================================

      module CAR_CU_Interface_Class
          implicit none

          private
          public :: construct_Interface_Class, destruct_Interface_Class

!===========================================================================
!
! CAR_CU_Interface Class type definition
! 
!===========================================================================

          type, public :: CAR_CU_Interface
              integer        :: this
              character * 60 :: file_name
              logical        :: verbose
              integer        :: rtt_format
          end type CAR_CU_Interface 

!===========================================================================
!
! Interfaces
! 
!===========================================================================

          interface construct_Interface_Class
              module procedure CAR_CU_Interface_construct
          end interface

          interface destruct_Interface_Class
              module procedure CAR_CU_Interface_destruct
          end interface

          contains

!===========================================================================
!
! Subroutines
! 
!===========================================================================
! Construct a new C++ CAR_CU_Interface class object (self).

              subroutine CAR_CU_Interface_construct(self)
                  type(CAR_CU_Interface), intent(inout) :: self
                  integer  :: bool_verbose = 0
                  if (self%verbose)  bool_verbose = 1

                  call construct_car_cu_interface(self%this, self%file_name,&
                      bool_verbose, self%rtt_format)

              end subroutine CAR_CU_Interface_construct

! Destroy a C++ CAR_CU_Interface class object (self).
              subroutine CAR_CU_Interface_destruct(self)
                  type(CAR_CU_Interface), intent(inout) :: self

                  call destruct_car_cu_interface(self%this)

              end subroutine CAR_CU_Interface_destruct

       end module CAR_CU_Interface_Class

!---------------------------------------------------------------------------
!                              end of amr_mesh/Shadow_CAR_CU_Interface.f90
!---------------------------------------------------------------------------
