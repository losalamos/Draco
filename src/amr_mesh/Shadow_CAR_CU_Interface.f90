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
          public :: construct_Interface, destruct_Interface,            &
                    get_surface_source_size, get_surface_source_pos,    &
                    get_surface_source_cells

!===========================================================================
! CAR_CU_Interface Class type definition
!===========================================================================

          type, public :: CAR_CU_Interface
              integer             :: this
              character*132       :: file_name
              logical             :: verbose
              integer             :: rtt_format
          end type CAR_CU_Interface 

!===========================================================================
! Interfaces
!===========================================================================

          interface construct_Interface
              module procedure CAR_CU_Interface_construct
          end interface

          interface destruct_Interface
              module procedure CAR_CU_Interface_destruct
          end interface

          interface get_surface_source_size
              module procedure get_ss_pos_size
              module procedure get_ss_cells_size
          end interface

          interface get_surface_source_pos
              module procedure get_ss_pos_set
          end interface

          interface get_surface_source_cells
              module procedure get_ss_cells_set
          end interface

          contains

!===========================================================================
! Constructors and destructors
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

!===========================================================================
! General interface accessor functions
!===========================================================================
! Return the number of sets of grouped surface source cells.
              function get_ss_pos_size(self)       result (data)
                  type(CAR_CU_Interface), intent(inout) :: self
                  integer                               :: data

                  call get_car_cu_ss_pos_size(self%this, data)

              end function get_ss_pos_size

! Return the number of grouped surface source cells in a given set.
              function get_ss_cells_size(self, surface)      result (data)
                  type(CAR_CU_Interface), intent(inout)   :: self
                  integer, intent(in)                     :: surface
                  integer                                 :: data

                  call get_car_cu_ss_cells_size(self%this, surface, data)

              end function get_ss_cells_size

! Return the position (lox, hix, etc.) of a given grouped surface source 
! cell set.
              function get_ss_pos_set(self, surface)       result (data)
                  type(CAR_CU_Interface), intent(inout) :: self
                  integer, intent(in)                   :: surface
                  character*3                           :: data

                  call get_car_cu_ss_cells_pos(self%this, surface, data)

              end function get_ss_pos_set

! Return the defined surface source cells in a given set.
              function get_ss_cells_set(self, surface, data_size) result (data)
                  type(CAR_CU_Interface), intent(inout) :: self
                  integer, intent(in)                   :: surface
                  integer, intent(in)                   :: data_size
                  integer, dimension(data_size)         :: data

                  call get_car_cu_ss_cells_set(self%this, surface, data, &
                                               data_size)

              end function get_ss_cells_set

!---------------------------------------------------------------------------
!                              end of CAR_CU_Interface module
!---------------------------------------------------------------------------

       end module CAR_CU_Interface_Class

!---------------------------------------------------------------------------
!                              end of amr_mesh/Shadow_CAR_CU_Interface.f90
!---------------------------------------------------------------------------
