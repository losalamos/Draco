!----------------------------------*-F90-*----------------------------------
! Shadow_CAR_CU_Mesh.f90
! B.T. Adams (bta@lanl.gov)
! 27 Sept 99
!---------------------------------------------------------------------------
! @> Shadow_CAR_CU_Mesh interface file
!---------------------------------------------------------------------------

!===========================================================================
! Shadow_CAR_CU_Mesh - 
!
! Purpose : Provides shadow interface functions for the C++ Continuous 
! Adaptive Refinement Cartesion Unstructured Mesh Class. Note that the class
! constructor is not shadowed because the mesh is constructed by the 
! Shadow_CAR_CU_Builder class object.
!
! revision history:
! -----------------
!  0) original
! 
!===========================================================================

      module CAR_CU_Mesh_Class
          implicit none

          private
          public :: destruct_Mesh_Class, get_dimension, get_num_cells,  &
                    get_num_nodes, get_num_corner_nodes,                &
                    get_num_face_nodes, get_num_adj, get_next_cell,     &
                    get_cell_volume, get_cell_face_area,                &
                    get_mesh_min_coord, get_mesh_max_coord,             &
                    get_cell_min_coord, get_cell_max_coord,             &
                    get_cell_mid_coord, get_cell_generation

!===========================================================================
!
! CAR_CU_Mesh Class type definition
! 
!===========================================================================

          type, public :: CAR_CU_Mesh
              integer             :: this
          end type CAR_CU_Mesh 

!===========================================================================
!
! Interfaces
! 
!===========================================================================

          interface destruct_Mesh_Class
              module procedure CAR_CU_Mesh_destruct
          end interface

          interface get_dimension
              module procedure get_dimension
          end interface

          interface get_num_cells
              module procedure get_num_cells
          end interface

          interface get_num_nodes
              module procedure get_num_nodes
          end interface

          interface get_num_corner_nodes
              module procedure get_num_corner_nodes
          end interface

          interface get_num_face_nodes
              module procedure get_num_face_nodes
          end interface
 
          interface get_num_adj
              module procedure get_num_adj
          end interface

          interface get_next_cell
              module procedure get_next_cell
          end interface

          interface get_cell_volume
              module procedure get_cell_volume
          end interface

          interface get_cell_face_area
              module procedure get_cell_face_area
          end interface

          interface get_cell_generation
              module procedure get_cell_generation
          end interface

          interface get_mesh_min_coord
              module procedure get_mesh_min_coord
          end interface

          interface get_mesh_max_coord
              module procedure get_mesh_max_coord
          end interface

          interface get_cell_min_coord
              module procedure get_cell_min_coord
          end interface

          interface get_cell_max_coord
              module procedure get_cell_max_coord
          end interface

          interface get_cell_mid_coord
              module procedure get_cell_mid_coord
          end interface

          contains

!===========================================================================
!
! Subroutines
! 
!===========================================================================
! Destroy a C++ CAR_CU_Mesh class object (self).
              subroutine CAR_CU_Mesh_destruct(self)
                  type(CAR_CU_Mesh), intent(inout) :: self

                  call destruct_car_cu_mesh(self%this)

              end subroutine CAR_CU_Mesh_destruct

!===========================================================================
!
! Functions
! 
!===========================================================================
! Return the geometry dimension of the mesh (self).
              integer function get_dimension(self) result(ndim)
                  type(CAR_CU_Mesh), intent(in) :: self
                  integer                       :: ndim

                  call get_mesh_dimension(self%this, ndim)

              end function get_dimension

! Return the number of cells in the mesh (self).
              integer function get_num_cells(self) result(ncells)
                  type(CAR_CU_Mesh), intent(in) :: self
                  integer                       :: ncells

                  call get_mesh_num_cells(self%this, ncells)

              end function get_num_cells


! Return the total number of nodes (nnodes) in the mesh (self).
              integer function get_num_nodes(self) result(nnodes)
                  type(CAR_CU_Mesh), intent(in) :: self
                  integer                       :: nnodes

                  call get_mesh_num_nodes(self%this, nnodes)

              end function get_num_nodes

! Return the number of cell-corner nodes (ncnodes) in the mesh (self).
              integer function get_num_corner_nodes(self) result(ncnodes)
                  type(CAR_CU_Mesh), intent(in) :: self
                  integer                       :: ncnodes

                  call get_mesh_num_corner_nodes(self%this,ncnodes)

              end function get_num_corner_nodes

! Return the number of face-centered nodes (nfnodes) in the mesh (self).
              integer function get_num_face_nodes(self) result(nfnodes)
                  type(CAR_CU_Mesh), intent(in) :: self
                  integer                       :: nfnodes

                  call get_mesh_num_face_nodes(self%this, nfnodes)

              end function get_num_face_nodes

! Return the number of cells that are adjacent to this cell face in the mesh
! (self).
              integer function get_num_adj(self, cell, face) result(num_adj)
                  type(CAR_CU_Mesh), intent(in) :: self
                  integer, intent(in)           :: cell, face
                  integer                       :: num_adj

                  call get_mesh_num_adj(self%this, cell, face, num_adj)

              end function get_num_adj

! Return the cell that is adjacent to this cell face in the mesh (self). An 
! optional index can be specified for cells that are adjacent to multiple 
! cells.
              integer function get_next_cell(self, cell, face, index)   &
                               result(adj_cell)
                  type(CAR_CU_Mesh), intent(in) :: self
                  integer, intent(in)           :: cell, face
                  integer                       :: adj_cell
                  integer,optional              :: index

                  if (.not. present(index)) then
                      call get_mesh_next_cell(self%this,cell,face,adj_cell)
                  else
                      call get_mesh_next_indexed_cell(self%this, cell,  &
                                                      face, index, adj_cell)
                  end if 

              end function get_next_cell

! Return the cell volume in the mesh (self).
              real function get_cell_volume(self, cell) result(volume)
                  type(CAR_CU_Mesh), intent(in) :: self
                  integer, intent(in)           :: cell
                  real                          :: volume

                  call get_mesh_cell_volume(self%this, cell, volume)

              end function get_cell_volume

! Return the cell face area in the mesh (self).
              real function get_cell_face_area(self, cell, face) result(area)
                  type(CAR_CU_Mesh), intent(in) :: self
                  integer, intent(in)           :: cell, face
                  real                          :: area

                  call get_mesh_cell_face_area(self%this, cell, face, area)

              end function get_cell_face_area

! Return the minimum coordinate value in a given direction for the mesh (self).
              real function get_mesh_min_coord(self, direction) result(min_val)
                  type(CAR_CU_Mesh), intent(in) :: self
                  integer, intent(in)           :: direction
                  real                          :: min_val

                  call get_mesh_min_coordinates(self%this, direction, min_val)

              end function get_mesh_min_coord

! Return the maximum coordinate value in a given direction for the mesh (self).
              real function get_mesh_max_coord(self, direction) result(max_val)
                  type(CAR_CU_Mesh), intent(in) :: self
                  integer, intent(in)           :: direction
                  real                          :: max_val

                  call get_mesh_max_coordinates(self%this, direction, max_val)

              end function get_mesh_max_coord

! Return the minimum coordinate value in a given direction for a cell in the 
! mesh (self).
              real function get_cell_min_coord(self, cell, dir) result(min_val)
                  type(CAR_CU_Mesh), intent(in) :: self
                  integer, intent(in)           :: cell, dir
                  real                          :: min_val

                  call get_mesh_cell_min_coord(self%this,cell,dir,min_val)

              end function get_cell_min_coord

! Return the maximum coordinate value in a given direction for a cell in the 
! mesh (self).
              real function get_cell_max_coord(self, cell, dir) result(max_val)
                  type(CAR_CU_Mesh), intent(in) :: self
                  integer, intent(in)           :: cell, dir
                  real                          :: max_val

                  call get_mesh_cell_max_coord(self%this,cell,dir,max_val)

              end function get_cell_max_coord

! Return the midpoint (i.e., center point) coordinate value in a given 
! direction for a cell in the mesh (self).
              real function get_cell_mid_coord(self, cell, dir) result(mid_val)
                  type(CAR_CU_Mesh), intent(in) :: self
                  integer, intent(in)           :: cell, dir
                  real                          :: mid_val

                  call get_mesh_cell_mid_coord(self%this,cell,dir,mid_val)

              end function get_cell_mid_coord

! Return the cell generation level in the mesh (self).
              integer function get_cell_generation(self, cell) result(gener)
                  type(CAR_CU_Mesh), intent(in) :: self
                  integer, intent(in)           :: cell
                  integer                       :: gener

                  call get_mesh_cell_generation(self%this, cell, gener)

              end function get_cell_generation

    end module CAR_CU_Mesh_Class


!---------------------------------------------------------------------------
!                              end of amr_mesh/Shadow_CAR_CU_Mesh.f90
!---------------------------------------------------------------------------
