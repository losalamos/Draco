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
!===========================================================================
! Constructors and destructors
!===========================================================================

          public :: destruct_Mesh_Class, construct_integer_CCSF_Class,  &
                    destruct_integer_CCSF_Class,                        &
                    construct_real_CCSF_Class, destruct_real_CCSF_Class

!===========================================================================
! General mesh scalar accessor functions
!===========================================================================

          public :: get_dimension, get_num_cells, get_num_nodes,        &
                    get_num_corner_nodes, get_num_face_nodes

!===========================================================================
! Layout accessor functions
!===========================================================================

          public :: get_num_adj, get_next_cell, get_cell_node,          &
                    get_cell_face_centered_node, get_cell_nodes,        &
                    get_cell_face_centered_nodes, get_cell_corner_nodes,&
                    get_cell_face_nodes

!===========================================================================
! Vertex accessor functions
!===========================================================================

          public :: get_vertices, get_corner_node_vertices,             &
                    get_face_centered_node_vertices, get_cell_vertices, &
                    get_cell_face_vertices, get_node_vertices

!===========================================================================
! Mesh geometry scalar accessor functions
!===========================================================================

          public :: get_cell_volume, get_cell_face_area,                &
                    get_mesh_min_coord, get_mesh_max_coord,             &
                    get_cell_min_coord, get_cell_mid_coord,             &
                    get_cell_max_coord, get_cell_width,                 &
                    get_cell_generation

!===========================================================================
! Mesh field accessor functions
!===========================================================================
! CCSF class objects

          public :: get_integer_CCSF, get_integer_CCSF_cell,            &
                    set_integer_CCSF, set_integer_CCSF_cell,            &
                    get_real_CCSF, get_real_CCSF_cell, set_real_CCSF,   &
                    set_real_CCSF_cell

!===========================================================================
! Class type definitions
!===========================================================================

          type, public :: CAR_CU_Mesh
              integer             :: this
          end type CAR_CU_Mesh 

          type, public :: integer_CCSF
              integer             :: this
              type(CAR_CU_Mesh)   :: mesh            
          end type integer_CCSF

          type, public :: real_CCSF
              integer             :: this
              type(CAR_CU_Mesh)   :: mesh
          end type real_CCSF

!===========================================================================
! Define interfaces
!===========================================================================

          interface destruct_Mesh_Class
              module procedure CAR_CU_Mesh_destruct
          end interface

          interface construct_integer_CCSF_Class
              module procedure CAR_CU_Mesh_int_CCSF_construct
          end interface

          interface destruct_integer_CCSF_Class
              module procedure CAR_CU_Mesh_int_CCSF_destruct
          end interface

          interface construct_real_CCSF_Class
              module procedure CAR_CU_Mesh_real_CCSF_construct
          end interface

          interface destruct_real_CCSF_Class
              module procedure CAR_CU_Mesh_real_CCSF_destruct
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

          interface get_cell_node
              module procedure get_cell_node
          end interface

          interface get_cell_face_centered_node
              module procedure get_cell_face_centered_node
          end interface

          interface get_cell_nodes
              module procedure get_cell_nodes
          end interface

          interface get_cell_face_centered_nodes
              module procedure get_cell_face_centered_nodes
          end interface

          interface get_cell_corner_nodes
              module procedure get_cell_corner_nodes
          end interface

          interface get_cell_face_nodes
              module procedure get_cell_face_nodes
          end interface

          interface get_vertices
              module procedure get_vertices
          end interface

          interface get_corner_node_vertices
              module procedure get_corner_node_vertices
          end interface

          interface get_face_centered_node_vertices
              module procedure get_face_centered_node_vertices
          end interface

          interface get_cell_vertices
              module procedure get_cell_vertices
          end interface

          interface get_cell_face_vertices
              module procedure get_cell_face_vertices
          end interface

          interface get_node_vertices
              module procedure get_node_vertices
          end interface

          interface get_cell_volume
              module procedure get_cell_volume
          end interface

          interface get_cell_face_area
              module procedure get_cell_face_area
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

          interface get_cell_mid_coord
              module procedure get_cell_mid_coord
          end interface

          interface get_cell_max_coord
              module procedure get_cell_max_coord
          end interface

          interface get_cell_width
              module procedure get_cell_width
          end interface

          interface get_cell_generation
              module procedure get_cell_generation
          end interface

          interface get_integer_CCSF
              module procedure get_integer_CCSF
          end interface

          interface get_integer_CCSF_cell
              module procedure get_integer_CCSF_cell
          end interface

          interface set_integer_CCSF
              module procedure set_integer_CCSF
          end interface

          interface set_integer_CCSF_cell
              module procedure set_integer_CCSF_cell
          end interface

          interface get_real_CCSF
              module procedure get_real_CCSF
          end interface

          interface get_real_CCSF_cell
              module procedure get_real_CCSF_cell
          end interface

          interface set_real_CCSF
              module procedure set_real_CCSF
          end interface

          interface set_real_CCSF_cell
              module procedure set_real_CCSF_cell
          end interface

          contains

!===========================================================================
! Constructors and destructors
!===========================================================================
! Destroy a C++ CAR_CU_Mesh class object (self).
              subroutine CAR_CU_Mesh_destruct(self)
                  type(CAR_CU_Mesh), intent(inout) :: self

                  call destruct_car_cu_mesh(self%this)

              end subroutine CAR_CU_Mesh_destruct

! Construct a C++ CAR_CU_Mesh integer CCSF class object (self). Initialization
! can be performed by including the optional data argument. An uninitialized
! CCSF is created if this argument is not specified.
              subroutine CAR_CU_Mesh_int_CCSF_construct(mesh, self, data)
                  type(CAR_CU_Mesh),  intent(in)           :: mesh
                  type(integer_CCSF), intent(inout)        :: self
                  integer, optional,                                    &
                           dimension(get_num_cells(mesh))  :: data
                  integer                                  :: data_size

                  if (.not. present(data)) then
                      call construct_mesh_ccsf_i(mesh%this,self%this)
                  else
                      data_size = get_num_cells(mesh)
                      call construct_mesh_ccsf_i_data(mesh%this,        &
                                     self%this, data, data_size)
                  endif
                  self%mesh = mesh

              end subroutine CAR_CU_Mesh_int_CCSF_construct

! Destroy a C++ CAR_CU_Mesh int CCSF class object (self).
              subroutine CAR_CU_Mesh_int_CCSF_destruct(self)
                  type(integer_CCSF), intent(inout) :: self

                  call destruct_mesh_ccsf_i(self%this)

              end subroutine CAR_CU_Mesh_int_CCSF_destruct

! Construct a C++ CAR_CU_Mesh real CCSF class object (self). Initialization
! can be performed by including the optional data argument. An uninitialized
! CCSF is created if this argument is not specified.
              subroutine CAR_CU_Mesh_real_CCSF_construct(mesh, self, data)
                  type(CAR_CU_Mesh), intent(in)            :: mesh
                  type(real_CCSF),   intent(inout)         :: self
                  real*8, optional,                                     &
                          dimension(get_num_cells(mesh))   :: data
                  integer                                  :: data_size

                  if (.not. present(data)) then
                      call construct_mesh_ccsf_d(mesh%this, self%this)
                  else
                      data_size = get_num_cells(mesh)
                      call construct_mesh_ccsf_d_data(mesh%this,        &
                                     self%this, data, data_size)
                  endif
                  self%mesh = mesh

              end subroutine CAR_CU_Mesh_real_CCSF_construct

! Destroy a C++ CAR_CU_Mesh real CCSF class object (self).
              subroutine CAR_CU_Mesh_real_CCSF_destruct(self)
                  type(real_CCSF), intent(inout) :: self

                  call destruct_mesh_ccsf_d(self%this)

              end subroutine CAR_CU_Mesh_real_CCSF_destruct

!===========================================================================
! General mesh scalar accessor functions
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

!===========================================================================
! Layout accessor functions
!===========================================================================
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
                      call get_mesh_next_specific_cell(self%this, cell,  &
                                                      face, index, adj_cell)
                  end if 

              end function get_next_cell

! Return the single cell node specified by the index.
              function get_cell_node(self, cell, node_index) result(node)
                  type(CAR_CU_Mesh), intent(in) :: self
                  integer, intent(in)           :: cell, node_index
                  integer                       :: node

                  call get_mesh_cell_node(self%this, cell, node_index, node)

              end function get_cell_node

! Return the single face-centered cell node specified by the face.
              function get_cell_face_centered_node(self, cell, face)     &
                  result(node)
                  type(CAR_CU_Mesh), intent(in) :: self
                  integer, intent(in)           :: cell, face
                  integer                       :: node

                  call get_mesh_cell_face_cen_node(self%this, cell, face, node)

              end function get_cell_face_centered_node


! Return an array of nodes that make up a cell, including both the corner nodes
! and the face-centered nodes.
              function get_cell_nodes(self, cell) result(nodes)
                  type(CAR_CU_Mesh), intent(in) :: self
                  integer, intent(in)           :: cell
                  integer, dimension(2**get_dimension(self) +            &
                                     2*get_dimension(self)) :: nodes
                  integer                       :: nsize

                  nsize = (2**get_dimension(self) + 2*get_dimension(self))
                  call get_mesh_cell_nodes(self%this, cell, nodes, nsize)

              end function get_cell_nodes

! Return an array of the cell face-centered nodes.
              function get_cell_face_centered_nodes(self, cell) result(nodes)
                  type(CAR_CU_Mesh), intent(in)              :: self
                  integer, intent(in)                        :: cell
                  integer, dimension(2*get_dimension(self))  :: nodes
                  integer                                    :: nsize

                  nsize = 2*get_dimension(self)
                  call get_mesh_cell_face_cen_nodes(self%this,cell,nodes,nsize)

              end function get_cell_face_centered_nodes

! Return an array of the cell corner nodes.
              function get_cell_corner_nodes(self, cell) result(nodes)
                  type(CAR_CU_Mesh), intent(in)              :: self
                  integer, intent(in)                        :: cell
                  integer, dimension(2**get_dimension(self)) :: nodes
                  integer                                    :: nsize

                  nsize = 2**get_dimension(self)
                  call get_mesh_cell_corner_nodes(self%this,cell,nodes,nsize)

              end function get_cell_corner_nodes

! Return an array of the nodes that comprise a cell face.
              function get_cell_face_nodes(self, cell, face) result(nodes)
                  type(CAR_CU_Mesh), intent(in)                   :: self
                  integer, intent(in)                             :: cell, face
                  integer, dimension(2**(get_dimension(self)-1))  :: nodes
                  integer                                         :: nsize

                  nsize = 2**(get_dimension(self) - 1)
                  call get_mesh_cell_face_nodes(self%this,cell,face,    &
                                                nodes,nsize)

              end function get_cell_face_nodes

!===========================================================================
! Vertex accessor functions
!===========================================================================
! Return the entire node vertex array (including both the corner and 
! face-centered nodes).

              function get_vertices(self) result(vertices)
                  type(CAR_CU_Mesh), intent(in)            :: self
                  real*8, dimension(get_num_nodes(self),                &
                                  get_dimension(self))     :: vertices
                  real*8, dimension(get_num_nodes(self) *               &
                                  get_dimension(self))     :: ret_vert
                  integer                                  :: node, dir, nsize

                  nsize = get_num_nodes(self) * get_dimension(self)
                  call get_mesh_vertices(self%this, ret_vert, nsize)

                  node = 1
                  do while (node .le. get_num_nodes(self))
                      dir = 1
                      do while (dir .le. get_dimension(self))
                          vertices(node, dir) =                         &
                              ret_vert(get_dimension(self)*(node-1) + dir)
                          dir = dir + 1
                      end do
                      node = node + 1
                  end do

              end function get_vertices

! Return an array containing the vertices for all of the cell corner nodes.

              function get_corner_node_vertices(self) result(vertices)
                  type(CAR_CU_Mesh), intent(in)             :: self
                  real*8, dimension(get_num_corner_nodes(self),         &
                                    get_dimension(self))    :: vertices
                  real*8, dimension(get_num_corner_nodes(self) *        &
                                    get_dimension(self))    :: ret_vert
                  integer                                   :: node, dir, nsize

                  nsize = get_num_corner_nodes(self) * get_dimension(self)
                  call get_mesh_corner_node_vertices(self%this,ret_vert,nsize)

                  node = 1
                  do while (node .le. get_num_corner_nodes(self))
                      dir = 1
                      do while (dir .le. get_dimension(self))
                          vertices(node, dir) =                         &
                              ret_vert(get_dimension(self)*(node-1) + dir)
                          dir = dir + 1
                      end do
                      node = node + 1
                  end do

              end function get_corner_node_vertices

! Return an array containing the vertices for all of the face_centered nodes.

              function get_face_centered_node_vertices(self) result(vertices)
                  type(CAR_CU_Mesh), intent(in)            :: self
                  real*8, dimension(get_num_face_nodes(self),           &
                                    get_dimension(self))   :: vertices
                  real*8, dimension(get_num_face_nodes(self) *          &
                                    get_dimension(self))   :: ret_vert
                  integer                                  :: node, dir, nsize

                  nsize = get_num_face_nodes(self) * get_dimension(self)
                  call get_mesh_face_cen_node_vertices(self%this,       &
                                                       ret_vert, nsize)

                  node = 1
                  do while (node .le. get_num_face_nodes(self))
                      dir = 1
                      do while (dir .le. get_dimension(self))
                          vertices(node, dir) =                         &
                              ret_vert(get_dimension(self)*(node-1) + dir)
                          dir = dir + 1
                      end do
                      node = node + 1
                  end do

              end function get_face_centered_node_vertices

! Return an array with all of a cell's vertices
              function get_cell_vertices(self, cell) result(vertices)
                  type(CAR_CU_Mesh), intent(in)             :: self
                  integer, intent(in)                       :: cell
                  real*8, dimension(2**get_dimension(self),             &
                                    get_dimension(self))    :: vertices
                  real*8, dimension(2**get_dimension(self) *            &
                                    get_dimension(self))    :: ret_vert

                  integer                                   :: node, nsize, dir

                  nsize = 2**get_dimension(self) * get_dimension(self)
                  call get_mesh_cell_vertices(self%this, cell,ret_vert, nsize)

                  node = 1
                  do while (node .le. 2**get_dimension(self))
                      dir = 1
                      do while (dir .le. get_dimension(self))
                          vertices(node, dir) =                         &
                              ret_vert(get_dimension(self)*(node-1) + dir)
                          dir = dir + 1
                      end do
                      node = node + 1
                  end do

              end function get_cell_vertices

! Return an array of all of the cell face vertices
              function get_cell_face_vertices(self,cell,face) result(vertices)
                  type(CAR_CU_Mesh), intent(in)             :: self
                  integer, intent(in)                       :: cell, face
                  real*8, dimension(2*(get_dimension(self)-1),          &
                                    get_dimension(self))    :: vertices
                  real*8, dimension(2*(get_dimension(self)-1) *         &
                                    get_dimension(self))    :: ret_vert
                  integer                                   :: node, nsize, dir

                  nsize = 2 * (get_dimension(self) - 1) * get_dimension(self)
                  call get_mesh_cell_face_vertices(self%this, cell,     &
                                                   face, ret_vert, nsize)
                  node = 1
                  do while (node .le. 2*(get_dimension(self) - 1))
                      dir = 1
                      do while (dir .le. get_dimension(self))
                          vertices(node, dir) =                         &
                              ret_vert(get_dimension(self)*(node-1) + dir)
                          dir = dir + 1
                      end do
                      node = node + 1
                  end do

              end function get_cell_face_vertices

! Return a single node's vertices
              function get_node_vertices(self, node) result(vertices)
                  type(CAR_CU_Mesh), intent(in)            :: self
                  integer, intent(in)                      :: node
                  real*8, dimension(get_dimension(self))   :: vertices
                  integer                                  :: nsize

                  nsize = get_dimension(self)
                  call get_mesh_node_vertices(self%this, node, vertices, nsize)

              end function get_node_vertices

!===========================================================================
! Mesh geometry scalar accessor functions
!===========================================================================
! Return the volume of the cell in the mesh (self).
                  function get_cell_volume(self, cell) result(volume)
                  type(CAR_CU_Mesh), intent(in) :: self
                  integer, intent(in)           :: cell
                  real*8                        :: volume

                  call get_mesh_cell_volume(self%this, cell, volume)

              end function get_cell_volume

! Return the face area of the cell in the mesh (self).
                  function get_cell_face_area(self, cell, face) result(area)
                  type(CAR_CU_Mesh), intent(in) :: self
                  integer, intent(in)           :: cell, face
                  real*8                        :: area

                  call get_mesh_cell_face_area(self%this, cell, face, area)

              end function get_cell_face_area

! Return the minimum coordinate value in a given direction for the mesh (self).
                  function get_mesh_min_coord(self, direction) result(min_val)
                  type(CAR_CU_Mesh), intent(in) :: self
                  integer, intent(in)           :: direction
                  real*8                        :: min_val

                  call get_mesh_min_coordinates(self%this, direction, min_val)

              end function get_mesh_min_coord

! Return the maximum coordinate value in a given direction for the mesh (self).
                  function get_mesh_max_coord(self, direction) result(max_val)
                  type(CAR_CU_Mesh), intent(in) :: self
                  integer, intent(in)           :: direction
                  real*8                        :: max_val

                  call get_mesh_max_coordinates(self%this, direction, max_val)

              end function get_mesh_max_coord

! Return the minimum coordinate value in a given direction for a cell in the 
! mesh (self).
                  function get_cell_min_coord(self, cell, dir) result(min_val)
                  type(CAR_CU_Mesh), intent(in) :: self
                  integer, intent(in)           :: cell, dir
                  real*8                        :: min_val

                  call get_mesh_cell_min_coord(self%this,cell,dir,min_val)

              end function get_cell_min_coord

! Return the midpoint (i.e., center point) coordinate value in a given 
! direction for a cell in the mesh (self).
                  function get_cell_mid_coord(self, cell, dir) result(mid_val)
                  type(CAR_CU_Mesh), intent(in) :: self
                  integer, intent(in)           :: cell, dir
                  real*8                        :: mid_val

                  call get_mesh_cell_mid_coord(self%this,cell,dir,mid_val)

              end function get_cell_mid_coord

! Return the maximum coordinate value in a given direction for a cell in the 
! mesh (self).
                  function get_cell_max_coord(self, cell, dir) result(max_val)
                  type(CAR_CU_Mesh), intent(in) :: self
                  integer, intent(in)           :: cell, dir
                  real*8                        :: max_val

                  call get_mesh_cell_max_coord(self%this,cell,dir,max_val)

              end function get_cell_max_coord

! Return the width in a given direction for a cell in the mesh (self).
                  function get_cell_width(self, cell, dir) result(width)
                  type(CAR_CU_Mesh), intent(in) :: self
                  integer, intent(in)           :: cell, dir
                  real*8                        :: width

                  call get_mesh_cell_width(self%this,cell,dir,width)

              end function get_cell_width

! Return the cell generation level in the mesh (self).
              integer function get_cell_generation(self, cell) result(gener)
                  type(CAR_CU_Mesh), intent(in) :: self
                  integer, intent(in)           :: cell
                  integer                       :: gener

                  call get_mesh_cell_generation(self%this, cell, gener)

              end function get_cell_generation
!===========================================================================
! Mesh field accessor functions
!===========================================================================
!===========================================================================
! integer CCSF class objects
!===========================================================================
! Return an entire C++ CAR_CU_Mesh integer CCSF class object (self).
              function get_integer_CCSF(self)    result(data)
                  type(integer_CCSF), intent(in)                :: self
                  integer, dimension(get_num_cells(self%mesh))  :: data
                  integer                                       :: data_size

                  data_size = get_num_cells(self%mesh)
                  call get_mesh_ccsf_i(self%mesh%this, self%this, data, &
                                       data_size)

              end function get_integer_CCSF

! Return a cell value from a C++ CAR_CU_Mesh integer CCSF class object (self).
              function get_integer_CCSF_cell(self, cell)      result(data)
                  type(integer_CCSF), intent(in)           :: self
                  integer, intent(in)                      :: cell
                  integer                                  :: data

                  call get_mesh_ccsf_i_cell(self%mesh%this, self%this,  &
                                            cell, data)

              end function get_integer_CCSF_cell

! Set an entire C++ CAR_CU_Mesh integer CCSF class object (self) (can also
! be done at initialization using the constructor).
              subroutine set_integer_CCSF(self, data)
                  type(integer_CCSF), intent(in)                :: self
                  integer, intent(in),                                  &
                           dimension(get_num_cells(self%mesh))  :: data
                  integer                                       :: data_size

                  data_size = get_num_cells(self%mesh)
                  call set_mesh_ccsf_i(self%mesh%this, self%this, data, &
                                       data_size)

              end subroutine set_integer_CCSF

! Set a cell value for a C++ CAR_CU_Mesh integer CCSF class object (self).
              subroutine set_integer_CCSF_cell(self, cell, data)
                  type(integer_CCSF), intent(in)           :: self
                  integer, intent(in)                      :: cell
                  integer, intent(in)                      :: data

                  call set_mesh_ccsf_i_cell(self%mesh%this, self%this,  &
                                            cell, data)

              end subroutine set_integer_CCSF_cell

!===========================================================================
! double CCSF class objects
!===========================================================================
! Return an entire C++ CAR_CU_Mesh double CCSF class object (self).
              function get_real_CCSF(self)       result(data)
                  type(real_CCSF), intent(in)                   :: self
                  real*8, dimension(get_num_cells(self%mesh))   :: data
                  integer                                       :: data_size

                  data_size = get_num_cells(self%mesh)
                  call get_mesh_ccsf_d(self%mesh%this, self%this, data, &
                                       data_size)

              end function get_real_CCSF

! Return a cell value from a C++ CAR_CU_Mesh double CCSF class object (self).
              function get_real_CCSF_cell(self, cell)     result(data)
                  type(real_CCSF), intent(in)             :: self
                  integer, intent(in)                     :: cell
                  real*8                                  :: data

                  call get_mesh_ccsf_d_cell(self%mesh%this, self%this,  &
                                            cell, data)

              end function get_real_CCSF_cell

! Set an entire C++ CAR_CU_Mesh double CCSF class object (self) (can also
! be done at initialization using the constructor).
              subroutine set_real_CCSF(self, data)
                  type(real_CCSF), intent(in)                   :: self
                  real*8, intent(in),                                   &
                           dimension(get_num_cells(self%mesh))  :: data
                  integer                                       :: data_size

                  data_size = get_num_cells(self%mesh)
                  call set_mesh_ccsf_d(self%mesh%this, self%this, data, &
                                       data_size)

              end subroutine set_real_CCSF

! Set a cell value from a C++ CAR_CU_Mesh double CCSF class object (self).
              subroutine set_real_CCSF_cell(self, cell, data)
                  type(real_CCSF), intent(in)              :: self
                  integer, intent(in)                      :: cell
                  real*8, intent(in)                       :: data

                  call set_mesh_ccsf_d_cell(self%mesh%this, self%this,  &
                                            cell, data)

              end subroutine set_real_CCSF_cell

    end module CAR_CU_Mesh_Class


!---------------------------------------------------------------------------
!                              end of amr_mesh/Shadow_CAR_CU_Mesh.f90
!---------------------------------------------------------------------------
