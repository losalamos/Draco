!----------------------------------*-F90-*----------------------------------
! Shadow_Builder.f90
! B.T. Adams (bta@lanl.gov)
! 27 Sept 99
!---------------------------------------------------------------------------
! @> Shadow_Builder interface file
!---------------------------------------------------------------------------

!===========================================================================
! Shadow_Builder - 
!
! Purpose : Provides shadow interface functions for the C++ Continuous 
! Adaptive Refinement Cartesion Unstructured Mesh Builder Class
!
! revision history:
! -----------------
!  0) original
! 
!===========================================================================

      module CAR_CU_Mesh_Builder_Class

          USE CAR_CU_Interface_Class

          implicit none

          private
          public :: construct_Mesh_Builder, destruct_Mesh_Builder

!===========================================================================
! CAR_CU_Mesh_Builder Class type definition
!===========================================================================

          type, public :: CAR_CU_Mesh_Builder
              integer        :: this
              integer        :: mesh
          end type CAR_CU_Mesh_Builder

!===========================================================================
! Interfaces
!===========================================================================

          interface construct_Mesh_Builder
              module procedure CAR_CU_Builder_construct
          end interface

          interface destruct_Mesh_Builder
              module procedure CAR_CU_Builder_destruct
          end interface

          contains

!===========================================================================
! Constructors and destructors
!===========================================================================
! Construct a new C++ CAR_CU_Mesh_Builder class object (self).

              subroutine CAR_CU_Builder_construct(self, itf_class)
                  type(CAR_CU_Mesh_Builder), intent(inout) :: self
                  type(CAR_CU_Interface), intent(inout)    :: itf_class
                  integer                                  :: bool_verbose = 0
                  if (itf_class%verbose)  bool_verbose = 1

                  call construct_car_cu_mesh_builder(self%this,         &
                       itf_class%this, itf_class%rtt_format,            &
                       bool_verbose, self%mesh)

              end subroutine CAR_CU_Builder_construct

! Destroy a C++ CAR_CU_Mesh_Builder class object (self).
              subroutine CAR_CU_Builder_destruct(self)
                  type(CAR_CU_Mesh_Builder), intent(inout) :: self

                  call destruct_car_cu_mesh_builder(self%this)

              end subroutine CAR_CU_Builder_destruct

!---------------------------------------------------------------------------
!                              end of CAR_CU_Mesh_Builder module
!---------------------------------------------------------------------------

       end module CAR_CU_Mesh_Builder_Class

!---------------------------------------------------------------------------
!                        end of amr_mesh_fort/Shadow_Builder.f90
!---------------------------------------------------------------------------
