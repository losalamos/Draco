!----------------------------------*-F90-*----------------------------------
! Shadow_Opacity.f90
! B.T. Adams (bta@lanl.gov)
! 27 Sept 99
!---------------------------------------------------------------------------
! @> Shadow_Opacity interface file
!---------------------------------------------------------------------------

!===========================================================================
! Shadow_Opacity - 
!
! Purpose : Provides shadow interface functions for the C++ Continuous 
! Adaptive Refinement Cartesion Unstructured Opacity Class. Note that the 
! class constructor is not shadowed because this object is constructed by 
! the CAR_CU_Opacity_Builder class object.
!
! revision history:
! -----------------
!  0) original
! 
!===========================================================================

      module CAR_CU_Opacity_Class
          implicit none

          private
!===========================================================================
! Constructors and destructors
!===========================================================================

          public ::  destruct_Opacity

!===========================================================================
! General Opacity scalar accessor functions
!===========================================================================

          public :: get_sigma_abs, get_sigma_thomson, get_planck,       &
                    get_fleck

!===========================================================================
! General Opacity scalar operator functions
!===========================================================================

          public :: get_fleck_planck, get_sigeffscat, get_sigeffabs

!===========================================================================
! Opacity field accessor functions
!===========================================================================

! none are defined.

!===========================================================================
! Class type definitions
!===========================================================================

          type, public :: CAR_CU_Opacity
              integer             :: this
              integer             :: mesh
          end type CAR_CU_Opacity 

!===========================================================================
! Define interfaces
!===========================================================================

          interface destruct_Opacity
              module procedure CAR_CU_Opacity_destruct
          end interface

          interface get_sigma_abs
              module procedure get_sigma_abs
          end interface

          interface get_sigma_thomson
              module procedure get_sigma_thomson
          end interface

          interface get_planck
              module procedure get_planck
          end interface

          interface get_fleck
              module procedure get_fleck
          end interface

          interface get_fleck_planck
              module procedure get_fleck_planck
          end interface

          interface get_sigeffscat
              module procedure get_sigeffscat
          end interface

          interface get_sigeffabs
              module procedure get_sigeffabs
          end interface

          contains

!===========================================================================
! Constructors and destructors
!===========================================================================
! Destroy a C++ CAR_CU_Opacity class object (self).
              subroutine CAR_CU_Opacity_destruct(self)
                  type(CAR_CU_Opacity), intent(inout) :: self

                  call destruct_car_cu_opacity(self%this)

              end subroutine CAR_CU_Opacity_destruct

!===========================================================================
! General Opacity scalar accessor functions
!===========================================================================
! Return sigma = kappa * rho for a cell in the opacity object (self).
              function get_sigma_abs(self, cell) result(data)
                  type(CAR_CU_Opacity), intent(in) :: self
                  integer, intent(in)              :: cell
                  real*8                           :: data

                  call get_car_cu_sigma_abs(self%this, self%mesh, cell, data)

              end function get_sigma_abs

! Return sigma_thomson = kappa_thomson * rho for a cell in the opacity object
! (self).
              function get_sigma_thomson(self, cell) result(data)
                  type(CAR_CU_Opacity), intent(in) :: self
                  integer, intent(in)              :: cell
                  real*8                           :: data

                  call get_car_cu_sigma_thomson(self%this,self%mesh,cell,data)

              end function get_sigma_thomson

! Return planckian opacity for a cell in the opacity object (self).
              function get_planck(self, cell) result(data)
                  type(CAR_CU_Opacity), intent(in) :: self
                  integer, intent(in)              :: cell
                  real*8                           :: data

                  call get_car_cu_planck(self%this, self%mesh, cell, data)

              end function get_planck

! Return fleck factor for a cell in the opacity object (self).
              function get_fleck(self, cell) result(data)
                  type(CAR_CU_Opacity), intent(in) :: self
                  integer, intent(in)              :: cell
                  real*8                           :: data

                  call get_car_cu_fleck(self%this, self%mesh, cell, data)

              end function get_fleck

!===========================================================================
! General Opacity scalar operator functions
!===========================================================================
! Return Fleck factor * Planckian opacity for a cell in the opacity object
! (self).
              function get_fleck_planck(self, cell) result(data)
                  type(CAR_CU_Opacity), intent(in) :: self
                  integer, intent(in)              :: cell
                  real*8                           :: data

                  call get_car_cu_fleck_planck(self%this,self%mesh,cell,data)

              end function get_fleck_planck

! Return the Fleck effective scatter cross section for a cell in the opacity 
! object (self).
              function get_sigeffscat(self, cell) result(data)
                  type(CAR_CU_Opacity), intent(in) :: self
                  integer, intent(in)              :: cell
                  real*8                           :: data

                  call get_car_cu_sigeffscat(self%this, self%mesh, cell, data)

              end function get_sigeffscat


! Return the Fleck effective absorption cross section for a cell in the 
! opacity object (self).
              function get_sigeffabs(self, cell) result(data)
                  type(CAR_CU_Opacity), intent(in) :: self
                  integer, intent(in)              :: cell
                  real*8                           :: data

                  call get_car_cu_sigeffabs(self%this, self%mesh, cell, data)

              end function get_sigeffabs

!---------------------------------------------------------------------------
!                              end of CAR_CU_Opacity_Class module
!---------------------------------------------------------------------------

      end module CAR_CU_Opacity_Class

!---------------------------------------------------------------------------
!                              end of amr_mesh_fort/Shadow_Opacity.f90
!---------------------------------------------------------------------------
