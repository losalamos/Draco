!----------------------------------*-F90-*----------------------------------
! Shadow_CAR_CU_Mesh_Builder.f90
! B.T. Adams (bta@lanl.gov)
! 27 Sept 99
!---------------------------------------------------------------------------
! @> Shadow_CAR_CU_Mesh_Builder interface file
!---------------------------------------------------------------------------

!===========================================================================
! Shadow_CAR_CU_Mesh_Builder - 
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
          public :: construct_Mesh_Builder, destruct_Mesh_Builder,      &
                    get_surface_source_size, get_surface_source_pos,    &
                    get_surface_source_cells

!===========================================================================
!
! CAR_CU_Mesh_Builder Class type definition
! 
!===========================================================================

          type, public :: CAR_CU_Mesh_Builder
              integer        :: this
              integer        :: mesh
          end type CAR_CU_Mesh_Builder

!===========================================================================
!
! Interfaces
! 
!===========================================================================

          interface construct_Mesh_Builder
              module procedure CAR_CU_Builder_construct
          end interface

          interface destruct_Mesh_Builder
              module procedure CAR_CU_Builder_destruct
          end interface

          interface get_surface_source_size
              module procedure get_ss_size
              module procedure get_ss_cells_size
          end interface

          interface get_surface_source_pos
              module procedure get_ss_set_pos
          end interface

          interface get_surface_source_cells
              module procedure get_ss_set_cells
          end interface

          contains

!===========================================================================
! Constructors and destructors
!===========================================================================
! Construct a new C++ CAR_CU_Mesh_Builder class object (self).

              subroutine CAR_CU_Builder_construct(self, itf_class)
                  type(CAR_CU_Mesh_Builder), intent(inout) :: self
                  type(CAR_CU_Interface), intent(inout)    :: itf_class
                  integer  :: bool_verbose = 0
                  if (itf_class%verbose)  bool_verbose = 1

                  call construct_car_cu_mesh_builder(self%this,         &
                      itf_class%this, itf_class%rtt_format,             &
                      bool_verbose, self%mesh)

              end subroutine CAR_CU_Builder_construct

! Destroy a C++ CAR_CU_Mesh_Builder class object (self).
              subroutine CAR_CU_Builder_destruct(self)
                  type(CAR_CU_Mesh_Builder), intent(inout) :: self

                  call destruct_car_cu_mesh_builder(self%this)

              end subroutine CAR_CU_Builder_destruct

!===========================================================================
! General mesh builder accessor functions
!===========================================================================
! Return the number of sets of grouped surface source cells.
              function get_ss_size(self)              result (data)
                  type(CAR_CU_Mesh_Builder), intent(inout) :: self
                  integer                                  :: data

                  call get_car_cu_ss_size(self%this, data)

              end function get_ss_size

! Return the number of grouped surface source cells in a given set.
              function get_ss_cells_size(self, surface) result (data)
                  type(CAR_CU_Mesh_Builder), intent(inout)   :: self
                  integer, intent(in)                        :: surface
                  integer                                    :: data

                  call get_car_cu_ss_set_size(self%this, surface, data)

              end function get_ss_cells_size

! Return the position (lox, hix, etc.) of a given grouped surface source 
! cell set.
              function get_ss_set_pos(self, surface)      result (data)
                  type(CAR_CU_Mesh_Builder), intent(inout) :: self
                  integer, intent(in)                      :: surface
                  character*3                              :: data

                  call get_car_cu_ss_set_pos(self%this, surface, data)

              end function get_ss_set_pos

! Return the defined surface source cells in a given set.
              function get_ss_set_cells(self, surface, data_size) result (data)
                  type(CAR_CU_Mesh_Builder), intent(inout) :: self
                  integer, intent(in)                      :: surface
                  integer, intent(in)                      :: data_size
                  integer, dimension(data_size)            :: data

                  call get_car_cu_ss_cell_set(self%this, surface, data, &
                                              data_size)

              end function get_ss_set_cells

!---------------------------------------------------------------------------
!                              end of CAR_CU_Mesh_Builder module
!---------------------------------------------------------------------------

       end module CAR_CU_Mesh_Builder_Class

!---------------------------------------------------------------------------
!                        end of amr_mesh/Shadow_CAR_CU_Mesh_Builder.f90
!---------------------------------------------------------------------------
