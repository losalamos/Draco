!----------------------------------*-F90-*----------------------------------
! Shadow_CAR_CU_Mat_State.f90
! B.T. Adams (bta@lanl.gov)
! 27 Sept 99
!---------------------------------------------------------------------------
! @> Shadow_CAR_CU_Mat_State interface file
!---------------------------------------------------------------------------

!===========================================================================
! Shadow_CAR_CU_Mat_State - 
!
! Purpose : Provides shadow interface functions for the C++ Continuous 
! Adaptive Refinement Cartesion Unstructured Mat_State Class. Note that the 
! class constructor is not shadowed because the mesh is constructed by the 
! Shadow_CAR_CU_Opacity_Builder class object.
!
! revision history:
! -----------------
!  0) original
! 
!===========================================================================

      module CAR_CU_Mat_State_Class
          implicit none

          private
!===========================================================================
! Constructors and destructors
!===========================================================================

          public ::  destruct_Mat_State

!===========================================================================
! General Mat_State scalar accessor functions
!===========================================================================

          public :: get_density, set_density, get_temperature,          &
                    set_temperature, get_dedt, set_dedt,                &
                    get_specific_heat

!===========================================================================
! Mat_State field accessor functions
!===========================================================================

! none are defined.

!===========================================================================
! Class type definitions
!===========================================================================

          type, public :: CAR_CU_Mat_State
              integer             :: this
              integer             :: mesh
          end type CAR_CU_Mat_State 

!===========================================================================
! Define interfaces
!===========================================================================

          interface destruct_Mat_State
              module procedure CAR_CU_Mat_State_destruct
          end interface

          interface get_density
              module procedure get_rho
          end interface

          interface set_density
              module procedure set_rho
          end interface

          interface get_temperature
              module procedure get_temp
          end interface

          interface set_temperature
              module procedure set_temp
          end interface

          interface get_dedt
              module procedure get_dedt
          end interface

          interface set_dedt
              module procedure set_dedt
          end interface

          interface get_specific_heat
              module procedure get_spec_heat
          end interface

          contains

!===========================================================================
! Constructors and destructors
!===========================================================================
! Destroy a C++ CAR_CU_Mat_State class object (self).
              subroutine CAR_CU_Mat_State_destruct(self)
                  type(CAR_CU_Mat_State), intent(inout) :: self

                  call destruct_car_cu_mat_state(self%this)

              end subroutine CAR_CU_Mat_State_destruct

!===========================================================================
! General Mat_State scalar accessor functions
!===========================================================================
! Return the cell density.
              function get_rho(self, cell)      result(data)
                  type(CAR_CU_Mat_State),intent(in) :: self
                  integer, intent(in)               :: cell
                  real*8                            :: data

                  call get_car_cu_rho(self%this, self%mesh, cell, data)

              end function get_rho

! Set the cell density.
              subroutine set_rho(self, cell, data)
                  type(CAR_CU_Mat_State),intent(in) :: self
                  integer, intent(in)               :: cell
                  real*8 , intent(in)               :: data

                  call set_car_cu_rho(self%this, self%mesh, cell, data)

              end subroutine set_rho

! Return the cell temperature.
              function get_temp(self, cell)      result(data)
                  type(CAR_CU_Mat_State),intent(in)  :: self
                  integer, intent(in)                :: cell
                  real*8                             :: data

                  call get_car_cu_temp(self%this, self%mesh, cell, data)

              end function get_temp

! Set the cell temperature.
              subroutine set_temp(self, cell, data)
                  type(CAR_CU_Mat_State),intent(in) :: self
                  integer, intent(in)               :: cell
                  real*8 , intent(in)               :: data

                  call set_car_cu_temp(self%this, self%mesh, cell, data)

              end subroutine set_temp

! Return the cell gradient of the cell energy with respect to temperature.
              function get_dedt(self, cell)      result(data)
                  type(CAR_CU_Mat_State),intent(in)  :: self
                  integer, intent(in)                :: cell
                  real*8                             :: data

                  call get_car_cu_dedt(self%this, self%mesh, cell, data)

              end function get_dedt

! Set the cell gradient of the cell energy with respect to temperature.
              subroutine set_dedt(self, cell, data)
                  type(CAR_CU_Mat_State),intent(in) :: self
                  integer, intent(in)               :: cell
                  real*8 , intent(in)               :: data

                  call set_car_cu_dedt(self%this, self%mesh, cell, data)

              end subroutine set_dedt

! Return the cell cell specific heat.
              function get_spec_heat(self, cell)      result(data)
                  type(CAR_CU_Mat_State),intent(in)       :: self
                  integer, intent(in)                     :: cell
                  real*8                                  :: data

                  call get_car_cu_spec_heat(self%this, self%mesh, cell, data)

              end function get_spec_heat

!---------------------------------------------------------------------------
!                              end of CAR_CU_Mat_State_Class module
!---------------------------------------------------------------------------

      end module CAR_CU_Mat_State_Class

!---------------------------------------------------------------------------
!                              end of amr_mesh/Shadow_CAR_CU_Mat_State.f90
!---------------------------------------------------------------------------
