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
                    get_surf_src_size, get_surf_src_pos,                &
                    get_surf_src_dist, get_surf_src_temperature,        &
                    get_surf_src_cells, get_vol_src, get_rad_src,       &
                    get_rad_src_tend, get_analytic_opacity,             &
                    get_analytic_specific_heat, get_implicitness,       &
                    get_time_step, get_capacity, get_num_cycles,        &
                    get_print_frequency

!===========================================================================
! CAR_CU_Interface Class type definition
!===========================================================================

          type, public :: CAR_CU_Interface
              integer                        :: this
              character*132                  :: file_name
              logical                        :: verbose
              integer                        :: rtt_format
              integer                        :: ss_size
              integer, dimension(:), pointer :: ss_cells
              integer                        :: ncells
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

          interface get_surf_src_size
              module procedure get_ss_pos_size
              module procedure get_ss_cells_size
          end interface

          interface get_surf_src_pos
              module procedure get_ss_pos
              module procedure get_ss_pos_set
          end interface

          interface get_surf_src_dist
              module procedure get_ss_dist
          end interface

          interface get_surf_src_temperature
              module procedure get_ss_temp
              module procedure get_ss_temp_set
          end interface

          interface get_surf_src_cells
              module procedure get_ss_cells_set
          end interface

          interface get_vol_src
              module procedure get_vol_src
              module procedure get_cell_vol_src
          end interface

          interface get_rad_src
              module procedure get_rad_src
              module procedure get_cell_rad_src
          end interface

          interface get_rad_src_tend
              module procedure get_rad_s_tend
          end interface

          interface get_analytic_opacity
              module procedure get_analy_opacity
          end interface

          interface get_analytic_specific_heat
              module procedure get_analy_spec_heat
          end interface

          interface get_implicitness
              module procedure get_implicit
          end interface

          interface get_time_step
              module procedure get_time_step
          end interface

          interface get_capacity
              module procedure get_capacity
          end interface

          interface get_num_cycles
              module procedure get_num_cycles
          end interface

          interface get_print_frequency
              module procedure get_print_freq
          end interface

          contains

!===========================================================================
! Constructors and destructors
!===========================================================================
! Construct a new C++ CAR_CU_Interface class object (self).

              subroutine CAR_CU_Interface_construct(self)
                  type(CAR_CU_Interface), intent(inout) :: self
                  integer                               :: bool_verbose = 0
                  integer :: test
                  integer                               :: surface
                  if (self%verbose)  bool_verbose = 1

                  call construct_car_cu_interface(self%this,            &
                           self%file_name, bool_verbose, self%rtt_format)

                  ! Set the number of surface sources
                  self%ss_size = get_surf_src_size(self)

                  ! Set the number of cells per surface source in the
                  ! interface class object
                  allocate(self%ss_cells(self%ss_size))

                  surface = 1
                  do while (surface .le. self%ss_size)
                      self%ss_cells = get_surf_src_size(self,surface)
                      surface = surface + 1
                  end do

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

! Return the position (lox, hix, etc.) of all of the grouped surface source 
! cells.
              function get_ss_pos(self)            result (data)
                  type(CAR_CU_Interface), intent(inout) :: self
                  character*3, dimension(self%ss_size)  :: data
                  character*(3 * self%ss_size)          :: ret_data
                  integer                               :: data_size
                  integer                               :: ss, index

                  data_size = 3 * self%ss_size
                  call get_car_cu_ss_pos(self%this, ret_data, data_size)

                  ss = 1
                  do while (ss .le. self%ss_size)
                      index = 1
                      do while (index .le. 3)
                          data(ss)(index:index) =                       &
                              ret_data(3*(ss -1) + index : 3*(ss -1) + index)
                          index = index + 1
                      end do
                      ss = ss + 1
                  end do

              end function get_ss_pos

! Return the position (lox, hix, etc.) of a given grouped surface source 
! cell set.
              function get_ss_pos_set(self, surface)       result (data)
                  type(CAR_CU_Interface), intent(inout) :: self
                  integer, intent(in)                   :: surface
                  character*3                           :: data

                  call get_car_cu_ss_cells_pos(self%this, surface, data)

              end function get_ss_pos_set

! Return the surface source angular distribution.
              function get_ss_dist(self)           result (data)
                  type(CAR_CU_Interface), intent(inout) :: self
                  character*6                           :: data

                  call get_car_cu_ss_dist(self%this, data)

              end function get_ss_dist

! Return the temperature of all of the grouped surface source cells.
              function get_ss_temp(self)       result (data)
                  type(CAR_CU_Interface), intent(inout) :: self
                  real*8, dimension(self%ss_size)       :: data
                  integer                               :: data_size

                  data_size = self%ss_size
                  call get_car_cu_ss_temp(self%this, data, data_size)

              end function get_ss_temp

! Return the temperature of a given grouped surface source cell set.
              function get_ss_temp_set(self, surface)       result (data)
                  type(CAR_CU_Interface), intent(inout) :: self
                  integer, intent(in)                   :: surface
                  real*8                                :: data

                  call get_car_cu_ss_cells_temp(self%this, surface, data)

              end function get_ss_temp_set

! Return the defined surface source cells in a given set.
              function get_ss_cells_set(self, surface)  result (data)
                  type(CAR_CU_Interface), intent(inout)      :: self
                  integer, intent(in)                        :: surface
                  integer, dimension(self%ss_cells(surface)) :: data
                  integer                                    :: data_size

                  data_size = self%ss_cells(surface)
                  call get_car_cu_ss_cells_set(self%this, surface, data, &
                                               data_size)

              end function get_ss_cells_set

! Return the mesh volumetric sources for all cells.
              function get_vol_src(self)           result (data)
                  type(CAR_CU_Interface), intent(inout) :: self
                  real*8, dimension(self%ncells)        :: data
                  integer                               :: data_size

                  data_size = self%ncells
                  call get_car_cu_vol_src(self%this, data, data_size)

              end function get_vol_src

! Return the volumetric source for a single cell.
              function get_cell_vol_src(self, cell) result (data)
                  type(CAR_CU_Interface), intent(inout)  :: self
                  integer, intent(in)                    :: cell
                  real*8                                 :: data

                  call get_car_cu_cell_vol_src(self%this, cell, data)

              end function get_cell_vol_src

! Return the mesh radiation sources for all cells.
              function get_rad_src(self)           result (data)
                  type(CAR_CU_Interface), intent(inout) :: self
                  real*8, dimension(self%ncells)        :: data
                  integer                               :: data_size

                  data_size = self%ncells
                  call get_car_cu_rad_src(self%this, data, data_size)

              end function get_rad_src

! Return the radiation source for a single cell.
              function get_cell_rad_src(self, cell) result (data)
                  type(CAR_CU_Interface), intent(inout)  :: self
                  integer, intent(in)                    :: cell
                  real*8                                 :: data

                  call get_car_cu_cell_rad_src(self%this, cell, data)

              end function get_cell_rad_src

! Return the cut-off time for a radiation source.
              function get_rad_s_tend(self)        result (data)
                  type(CAR_CU_Interface), intent(inout) :: self
                  real*8                                :: data

                  call get_car_cu_rad_s_tend(self%this, data)

              end function get_rad_s_tend

! Return the mesh analytic opacity type
              function get_analy_opacity(self) result(data)
                  type(CAR_CU_Interface),intent(in)  :: self
                  character*8                        :: data

                  call get_car_cu_analy_opacity(self%this, data)

              end function get_analy_opacity

! Return the mesh analytic specific heat type
              function get_analy_spec_heat(self) result(data)
                  type(CAR_CU_Interface),intent(in)  :: self
                  character*8                        :: data

                  call get_car_cu_analy_sp_heat(self%this, data)

              end function get_analy_spec_heat

! Return the mesh implicitness factor (Fleck's alpha).
              function get_implicit(self)        result(data)
                  type(CAR_CU_Interface),intent(in)  :: self
                  integer                            :: data

                  call get_car_cu_implicit(self%this, data)

              end function get_implicit

! Return the initial time step size.
              function get_time_step(self)       result(data)
                  type(CAR_CU_Interface),intent(in)  :: self
                  real*8                             :: data

                  call get_car_cu_time_step(self%this, data)

              end function get_time_step

! Return the processor capacity (cells/processor).
              function get_capacity(self)        result(data)
                  type(CAR_CU_Interface),intent(in)  :: self
                  integer                            :: data

                  call get_car_cu_capacity(self%this, data)

              end function get_capacity

! Return the number of cycles to run.
              function get_num_cycles(self)        result(data)
                  type(CAR_CU_Interface),intent(in)  :: self
                  integer                            :: data

                  call get_car_cu_num_cycles(self%this, data)

              end function get_num_cycles

! Return the printout frequency.
              function get_print_freq(self)        result(data)
                  type(CAR_CU_Interface),intent(in)  :: self
                  integer                            :: data

                  call get_car_cu_print_freq(self%this, data)

              end function get_print_freq

!---------------------------------------------------------------------------
!                              end of CAR_CU_Interface module
!---------------------------------------------------------------------------

       end module CAR_CU_Interface_Class

!---------------------------------------------------------------------------
!                              end of amr_mesh/Shadow_CAR_CU_Interface.f90
!---------------------------------------------------------------------------
