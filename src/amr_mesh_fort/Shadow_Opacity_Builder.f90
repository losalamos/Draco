!----------------------------------*-F90-*----------------------------------
! Shadow_Opacity_Builder.f90
! B.T. Adams (bta@lanl.gov)
! 27 Sept 99
!---------------------------------------------------------------------------
! @> Shadow_Opacity_Builder interface file
!---------------------------------------------------------------------------

!===========================================================================
! Shadow_Opacity_Builder - 
!
! Purpose : Provides shadow interface functions for the C++ Continuous 
! Adaptive Refinement Cartesion Unstructured Opacity Builder Class
!
! revision history:
! -----------------
!  0) original
! 
!===========================================================================

      module CAR_CU_Opacity_Builder_Class

          USE CAR_CU_Interface_Class
          USE CAR_CU_Mesh_Class

          implicit none

          private
          public :: construct_Opacity_Builder, destruct_Opacity_Builder

!===========================================================================
!
! Opacity_Builder Class type definition
! 
!===========================================================================

          type, public :: CAR_CU_Opacity_Builder
              integer        :: this
              integer        :: mesh
              integer        :: matl_state
              integer        :: opacity
          end type CAR_CU_Opacity_Builder 

!===========================================================================
!
! Interfaces
! 
!===========================================================================

          interface construct_Opacity_Builder
              module procedure Opacity_Builder_construct
          end interface

          interface destruct_Opacity_Builder
              module procedure Opacity_Builder_destruct
          end interface

          contains

!===========================================================================
!
! Subroutines
! 
!===========================================================================
! Construct a new C++ Opacity_Builder class object (self).

              subroutine Opacity_Builder_construct(self, itf_class, mesh_class)
                  type(CAR_CU_Opacity_Builder), intent(inout) :: self
                  type(CAR_CU_Interface), intent(in)          :: itf_class
                  type(CAR_CU_Mesh), intent(in)               :: mesh_class
                  integer  :: bool_verbose = 0
                  if (itf_class%verbose)  bool_verbose = 1

                  call construct_car_cu_opac_builder(self%this,           &
                      itf_class%this, mesh_class%this, bool_verbose,      &
                      self%matl_state, self%opacity)

                  self%mesh = mesh_class%this

              end subroutine Opacity_Builder_construct

! Destroy a C++ Opacity_Builder class object (self).
              subroutine Opacity_Builder_destruct(self)
                  type(CAR_CU_Opacity_Builder), intent(inout) :: self

                  call destruct_car_cu_opac_builder(self%this)

              end subroutine Opacity_Builder_destruct

       end module CAR_CU_Opacity_Builder_Class

!---------------------------------------------------------------------------
!                       end of amr_mesh/Shadow_CAR_CU_Opacity_Builder.f90
!---------------------------------------------------------------------------
