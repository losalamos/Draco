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

          public :: destruct_Mesh_Class, construct_CCSF_Class,          &
                    destruct_CCSF_Class, construct_CCVF_Class,          &
                    destruct_CCVF_Class, construct_FCSF_Class,          &
                    destruct_FCSF_Class, construct_FCDSF_Class,         &
                    destruct_FCDSF_Class

!===========================================================================
! General mesh scalar accessor functions
!===========================================================================

          public :: get_num_dims, get_num_cells, get_num_nodes,         &
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

          public :: get_CCSF, set_CCSF, get_CCVF, set_CCVF, get_FCSF,   &
                    set_FCSF, get_FCDSF, set_FCDSF

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

          type, public :: integer_FCSF
              integer             :: this
              type(CAR_CU_Mesh)   :: mesh            
          end type integer_FCSF

          type, public :: real_FCSF
              integer             :: this
              type(CAR_CU_Mesh)   :: mesh
          end type real_FCSF

          type, public :: integer_FCDSF
              integer             :: this
              type(CAR_CU_Mesh)   :: mesh            
          end type integer_FCDSF

          type, public :: real_FCDSF
              integer             :: this
              type(CAR_CU_Mesh)   :: mesh
          end type real_FCDSF

          type, public :: integer_CCVF
              integer             :: this
              type(CAR_CU_Mesh)   :: mesh            
              integer             :: vec_size
          end type integer_CCVF

          type, public :: real_CCVF
              integer             :: this
              type(CAR_CU_Mesh)   :: mesh
              integer             :: vec_size
          end type real_CCVF


!===========================================================================
! Define interfaces
!===========================================================================

          interface destruct_Mesh_Class
              module procedure CAR_CU_Mesh_destruct
          end interface

          interface construct_CCSF_Class
              module procedure int_CCSF_construct
              module procedure real_CCSF_construct
          end interface

          interface destruct_CCSF_Class
              module procedure int_CCSF_destruct
              module procedure real_CCSF_destruct
          end interface

          interface construct_CCVF_Class
              module procedure int_CCVF_construct
              module procedure int_CCVF_arb_construct
              module procedure real_CCVF_construct
              module procedure real_CCVF_arb_construct
          end interface

          interface destruct_CCVF_Class
              module procedure int_CCVF_destruct
              module procedure real_CCVF_destruct
          end interface

          interface construct_FCSF_Class
              module procedure int_FCSF_construct
              module procedure real_FCSF_construct
          end interface

          interface destruct_FCSF_Class
              module procedure int_FCSF_destruct
              module procedure real_FCSF_destruct
          end interface

          interface construct_FCDSF_Class
              module procedure int_FCDSF_construct
              module procedure real_FCDSF_construct
          end interface

          interface destruct_FCDSF_Class
              module procedure int_FCDSF_destruct
              module procedure real_FCDSF_destruct
          end interface

          interface get_num_dims
              module procedure get_num_dims
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

          interface get_CCSF
              module procedure get_integer_CCSF_all
              module procedure get_integer_CCSF_cell
              module procedure get_real_CCSF_all
              module procedure get_real_CCSF_cell
          end interface

          interface set_CCSF
              module procedure set_integer_CCSF_all
              module procedure set_integer_CCSF_cell
              module procedure set_real_CCSF_all
              module procedure set_real_CCSF_cell
          end interface

          interface get_CCVF
              module procedure get_integer_CCVF_all
              module procedure get_integer_CCVF_cell
              module procedure get_integer_CCVF_cell_dim
              module procedure get_real_CCVF_all
              module procedure get_real_CCVF_cell
              module procedure get_real_CCVF_cell_dim
          end interface

          interface set_CCVF
              module procedure set_integer_CCVF_all
              module procedure set_integer_CCVF_cell
              module procedure set_integer_CCVF_cell_dim
              module procedure set_real_CCVF_all
              module procedure set_real_CCVF_cell
              module procedure set_real_CCVF_cell_dim
          end interface

          interface get_FCSF
              module procedure get_integer_FCSF_all
              module procedure get_integer_FCSF_cell
              module procedure get_integer_FCSF_cell_face
              module procedure get_real_FCSF_all
              module procedure get_real_FCSF_cell
              module procedure get_real_FCSF_cell_face
          end interface

          interface set_FCSF
              module procedure set_integer_FCSF_all
              module procedure set_integer_FCSF_cell
              module procedure set_integer_FCSF_cell_face
              module procedure set_real_FCSF_all
              module procedure set_real_FCSF_cell
              module procedure set_real_FCSF_cell_face
          end interface

          interface get_FCDSF
              module procedure get_integer_FCDSF_all
              module procedure get_integer_FCDSF_cell
              module procedure get_integer_FCDSF_cell_face
              module procedure get_real_FCDSF_all
              module procedure get_real_FCDSF_cell
              module procedure get_real_FCDSF_cell_face
          end interface

          interface set_FCDSF
              module procedure set_integer_FCDSF_all
              module procedure set_integer_FCDSF_cell
              module procedure set_integer_FCDSF_cell_face
              module procedure set_real_FCDSF_all
              module procedure set_real_FCDSF_cell
              module procedure set_real_FCDSF_cell_face
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
              subroutine int_CCSF_construct(mesh, self, data)
                  type(CAR_CU_Mesh),  intent(in)           :: mesh
                  type(integer_CCSF), intent(inout)        :: self
                  integer, intent(in), optional,                        &
                           dimension(get_num_cells(mesh))  :: data
                  integer                                  :: data_size = 0

                  if (present(data)) data_size = get_num_cells(mesh)

                  call construct_mesh_ccsf_i(mesh%this, self%this,      &
                                             data, data_size)
                  self%mesh = mesh

              end subroutine int_CCSF_construct

! Destroy a C++ CAR_CU_Mesh int CCSF class object (self).
              subroutine int_CCSF_destruct(self)
                  type(integer_CCSF), intent(inout) :: self

                  call destruct_mesh_ccsf_i(self%this)

              end subroutine int_CCSF_destruct

! Construct a C++ CAR_CU_Mesh real CCSF class object (self). Initialization
! can be performed by including the optional data argument. An uninitialized
! CCSF is created if this argument is not specified.
              subroutine real_CCSF_construct(mesh, self, data)
                  type(CAR_CU_Mesh), intent(in)            :: mesh
                  type(real_CCSF),   intent(inout)         :: self
                  real*8, intent(in), optional,                         &
                          dimension(get_num_cells(mesh))   :: data
                  integer                                  :: data_size = 0

                  if (present(data)) data_size = get_num_cells(mesh)
                  call construct_mesh_ccsf_d(mesh%this, self%this,      &
                                                 data, data_size)
                  self%mesh = mesh

              end subroutine real_CCSF_construct

! Destroy a C++ CAR_CU_Mesh real CCSF class object (self).
              subroutine real_CCSF_destruct(self)
                  type(real_CCSF), intent(inout) :: self

                  call destruct_mesh_ccsf_d(self%this)

              end subroutine real_CCSF_destruct

! Construct a C++ CAR_CU_Mesh integer CCVF class object (self). Initialization
! can be performed by including the optional data argument. An uninitialized
! CCVF is created if this argument is not specified. This constructor defaults
! to the dimension of the second array element being the same as that of the 
! problem geometry.
              subroutine int_CCVF_construct(mesh, self, data)
                  type(CAR_CU_Mesh),  intent(in)           :: mesh
                  type(integer_CCVF), intent(inout)        :: self
                  integer, intent(in), optional,                        &
                           dimension(get_num_cells(mesh),               &
                                     get_num_dims(mesh))   :: data
                  integer, dimension(get_num_cells(mesh) *              &
                                     get_num_dims(mesh))   :: ret_data
                  integer                                  :: data_size = 0
                  integer                                  :: cell, dir
                  integer                                  :: vec_size

                  if (present(data)) then
                      cell = 1
                      do while(cell .le. get_num_cells(mesh))
                          dir = 1
                          do while (dir .le. get_num_dims(mesh))
                              ret_data(get_num_dims(mesh) * (cell-1) +  &
                                       dir) = data(cell, dir)
                              dir = dir + 1
                          end do
                          cell = cell + 1
                      end do
                      data_size = get_num_cells(mesh) * get_num_dims(mesh)
                  end if

                  vec_size = get_num_dims(mesh)
                  call construct_mesh_ccvf_i(mesh%this, self%this,      &
                                             ret_data, data_size, vec_size)
                  self%mesh = mesh
                  self%vec_size = vec_size

              end subroutine int_CCVF_construct

! Construct a C++ CAR_CU_Mesh integer CCVF class object (self). Initialization
! can be performed by including the optional data argument. An uninitialized
! CCVF is created if this argument is not specified. This constructor allows
! specification for the size of the second array index size.

              subroutine int_CCVF_arb_construct(mesh, self, vec_size, data)
                  type(CAR_CU_Mesh),  intent(in)           :: mesh
                  type(integer_CCVF), intent(inout)        :: self
                  integer,  intent(in)                     :: vec_size
                  integer, intent(in), optional,                        &
                           dimension(get_num_cells(mesh),vec_size) :: data
                  integer, dimension(get_num_cells(mesh)*vec_size) :: ret_data
                  integer                                  :: data_size = 0
                  integer                                  :: cell, dir

                  if (present(data)) then
                      cell = 1
                      do while(cell .le. get_num_cells(mesh))
                          dir = 1
                          do while (dir .le. get_num_dims(mesh))
                              ret_data(get_num_dims(mesh) * (cell-1) +  &
                                       dir) = data(cell, dir)
                              dir = dir + 1
                          end do
                          cell = cell + 1
                      end do
                      data_size = get_num_cells(mesh) * vec_size
                  end if

                  call construct_mesh_ccvf_i(mesh%this, self%this,      &
                                             ret_data, data_size, vec_size)
                  self%mesh = mesh
                  self%vec_size = vec_size

              end subroutine int_CCVF_arb_construct

! Destroy a C++ CAR_CU_Mesh int CCVF class object (self).
              subroutine int_CCVF_destruct(self)
                  type(integer_CCVF), intent(inout) :: self

                  call destruct_mesh_ccvf_i(self%this)

              end subroutine int_CCVF_destruct

! Construct a C++ CAR_CU_Mesh real CCVF class object (self). Initialization
! can be performed by including the optional data argument. An uninitialized
! CCVF is created if this argument is not specified. This constructor defaults
! to the dimension of the second array element being the same as that of the 
! problem geometry.
              subroutine real_CCVF_construct(mesh, self, data)
                  type(CAR_CU_Mesh), intent(in)            :: mesh
                  type(real_CCVF),   intent(inout)         :: self
                  real*8, intent(in), optional,                         &
                          dimension(get_num_cells(mesh),                &
                                    get_num_dims(mesh))    :: data
                  real*8, dimension(get_num_cells(mesh) *               &
                                    get_num_dims(mesh))    :: ret_data
                  integer                                  :: data_size = 0
                  integer                                  :: cell, dir
                  integer                                  :: vec_size

                  if (present(data)) then
                      cell = 1
                      do while(cell .le. get_num_cells(mesh))
                          dir = 1
                          do while (dir .le. get_num_dims(mesh))
                              ret_data(get_num_dims(mesh) * (cell-1) +  &
                                       dir) = data(cell, dir)
                              dir = dir + 1
                          end do
                          cell = cell + 1
                      end do

                      data_size = get_num_cells(mesh) * get_num_dims(mesh)
                  end if

                  vec_size = get_num_dims(mesh)
                  call construct_mesh_ccvf_d(mesh%this, self%this,      &
                                             ret_data, data_size, vec_size)
                  self%mesh = mesh
                  self%vec_size = vec_size

              end subroutine real_CCVF_construct

! Construct a C++ CAR_CU_Mesh real CCVF class object (self). Initialization
! can be performed by including the optional data argument. An uninitialized
! CCVF is created if this argument is not specified. This constructor allows
! specification for the size of the second array index size.
              subroutine real_CCVF_arb_construct(mesh, self, vec_size, data)
                  type(CAR_CU_Mesh), intent(in)            :: mesh
                  type(real_CCVF),   intent(inout)         :: self
                  integer, intent(in)                      :: vec_size
                  real*8, intent(in), optional,                         &
                          dimension(get_num_cells(mesh), vec_size)  :: data
                  real*8, dimension(get_num_cells(mesh) * vec_size) :: ret_data
                  integer                                  :: data_size = 0
                  integer                                  :: cell, dir

                  if (present(data)) then
                      cell = 1
                      do while(cell .le. get_num_cells(mesh))
                          dir = 1
                          do while (dir .le. get_num_dims(mesh))
                              ret_data(get_num_dims(mesh) * (cell-1) +  &
                                       dir) = data(cell, dir)
                              dir = dir + 1
                          end do
                          cell = cell + 1
                      end do

                      data_size = get_num_cells(mesh) * vec_size
                  end if

                  call construct_mesh_ccvf_d(mesh%this, self%this,      &
                                             ret_data, data_size, vec_size)
                  self%mesh = mesh
                  self%vec_size = vec_size

              end subroutine real_CCVF_arb_construct

! Destroy a C++ CAR_CU_Mesh real CCVF class object (self).
              subroutine real_CCVF_destruct(self)
                  type(real_CCVF), intent(inout) :: self

                  call destruct_mesh_ccvf_d(self%this)

              end subroutine real_CCVF_destruct

! Construct a C++ CAR_CU_Mesh integer FCSF class object (self). Initialization
! can be performed by including the optional data argument. An uninitialized
! FCSF is created if this argument is not specified.
              subroutine int_FCSF_construct(mesh, self, data)
                  type(CAR_CU_Mesh),  intent(in)               :: mesh
                  type(integer_FCSF), intent(inout)            :: self
                  integer, optional,                                    &
                           dimension(get_num_face_nodes(mesh)) :: data
                  integer                                      :: data_size = 0

                  if (present(data)) data_size = get_num_face_nodes(mesh)
                  call construct_mesh_fcsf_i(mesh%this, self%this,      &
                                             data, data_size)
                  self%mesh = mesh

              end subroutine int_FCSF_construct

! Destroy a C++ CAR_CU_Mesh int FCSF class object (self).
              subroutine int_FCSF_destruct(self)
                  type(integer_FCSF), intent(inout) :: self

                  call destruct_mesh_fcsf_i(self%this)

              end subroutine int_FCSF_destruct

! Construct a C++ CAR_CU_Mesh real FCSF class object (self). Initialization
! can be performed by including the optional data argument. An uninitialized
! FCSF is created if this argument is not specified.
              subroutine real_FCSF_construct(mesh, self, data)
                  type(CAR_CU_Mesh), intent(in)               :: mesh
                  type(real_FCSF),   intent(inout)            :: self
                  real*8, optional,                                     &
                          dimension(get_num_face_nodes(mesh)) :: data
                  integer                                     :: data_size = 0
                  integer                                     :: cell,face

                  if (present(data)) data_size = get_num_face_nodes(mesh)

                  call construct_mesh_fcsf_d(mesh%this, self%this,      &
                                             data, data_size)
                  self%mesh = mesh

              end subroutine real_FCSF_construct

! Destroy a C++ CAR_CU_Mesh real FCSF class object (self).
              subroutine real_FCSF_destruct(self)
                  type(real_FCSF), intent(inout) :: self

                  call destruct_mesh_fcsf_d(self%this)

              end subroutine real_FCSF_destruct

! Construct a C++ CAR_CU_Mesh integer FCDSF class object (self). Initialization
! can be performed by including the optional data argument. An uninitialized
! FCDSF is created if this argument is not specified.
              subroutine int_FCDSF_construct(mesh, self, data)
                  type(CAR_CU_Mesh),  intent(in)               :: mesh
                  type(integer_FCDSF), intent(inout)           :: self
                  integer, intent(in), optional,                        &
                           dimension(get_num_cells(mesh),               &
                                     2 * get_num_dims(mesh))   :: data
                  integer, dimension(get_num_cells(mesh) *              &
                                     2 * get_num_dims(mesh))   :: ret_data
                  integer                                      :: data_size = 0
                  integer                                      :: cell, face

                  if (present(data)) then
                      cell = 1
                      do while(cell .le. get_num_cells(mesh))
                          face = 1
                          do while (face .le. 2 * get_num_dims(mesh))
                              ret_data(2 * get_num_dims(mesh) *         &
                                  (cell-1) + face) = data(cell, face)
                              face = face + 1
                          end do
                          cell = cell + 1
                      end do

                      data_size = get_num_cells(mesh) * 2 * get_num_dims(mesh)
                  end if

                  call construct_mesh_fcdsf_i(mesh%this, self%this,     &
                                              ret_data, data_size)
                  self%mesh = mesh

              end subroutine int_FCDSF_construct

! Destroy a C++ CAR_CU_Mesh int FCDSF class object (self).
              subroutine int_FCDSF_destruct(self)
                  type(integer_FCDSF), intent(inout) :: self

                  call destruct_mesh_fcdsf_i(self%this)

              end subroutine int_FCDSF_destruct

! Construct a C++ CAR_CU_Mesh real FCDSF class object (self). Initialization
! can be performed by including the optional data argument. An uninitialized
! FCDSF is created if this argument is not specified.
              subroutine real_FCDSF_construct(mesh, self, data)
                  type(CAR_CU_Mesh), intent(in)               :: mesh
                  type(real_FCDSF),   intent(inout)           :: self
                  real*8, intent(in), optional,                         &
                          dimension(get_num_cells(mesh),                &
                                    2 * get_num_dims(mesh))   :: data
                  real*8, dimension(get_num_cells(mesh) *               &
                                    2 * get_num_dims(mesh))   :: ret_data
                  integer                                     :: data_size = 0
                  integer                                     :: cell, face

                  if (present(data)) then
                      cell = 1
                      do while(cell .le. get_num_cells(mesh))
                          face = 1
                          do while (face .le. 2 * get_num_dims(mesh))
                              ret_data(2 * get_num_dims(mesh) *         &
                                  (cell-1) + face) = data(cell, face)
                              face = face + 1
                          end do
                          cell = cell + 1
                      end do

                      data_size = get_num_cells(mesh) * 2 * get_num_dims(mesh)
                  end if

                  call construct_mesh_fcdsf_d(mesh%this, self%this,     &
                                              ret_data, data_size)
                  self%mesh = mesh

              end subroutine real_FCDSF_construct

! Destroy a C++ CAR_CU_Mesh real FCDSF class object (self).
              subroutine real_FCDSF_destruct(self)
                  type(real_FCDSF), intent(inout) :: self

                  call destruct_mesh_fcdsf_d(self%this)

              end subroutine real_FCDSF_destruct

!===========================================================================
! General mesh scalar accessor functions
!===========================================================================
! Return the geometry dimension of the mesh (self).
              integer function get_num_dims(self) result(ndim)
                  type(CAR_CU_Mesh), intent(in) :: self
                  integer                       :: ndim

                  call get_mesh_dimension(self%this, ndim)

              end function get_num_dims

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
                  integer, dimension(2**get_num_dims(self) +             &
                                     2*get_num_dims(self)) :: nodes
                  integer                       :: nsize

                  nsize = (2**get_num_dims(self) + 2*get_num_dims(self))
                  call get_mesh_cell_nodes(self%this, cell, nodes, nsize)

              end function get_cell_nodes

! Return an array of the cell face-centered nodes.
              function get_cell_face_centered_nodes(self, cell) result(nodes)
                  type(CAR_CU_Mesh), intent(in)              :: self
                  integer, intent(in)                        :: cell
                  integer, dimension(2*get_num_dims(self))   :: nodes
                  integer                                    :: nsize

                  nsize = 2*get_num_dims(self)
                  call get_mesh_cell_face_cen_nodes(self%this,cell,nodes,nsize)

              end function get_cell_face_centered_nodes

! Return an array of the cell corner nodes.
              function get_cell_corner_nodes(self, cell) result(nodes)
                  type(CAR_CU_Mesh), intent(in)              :: self
                  integer, intent(in)                        :: cell
                  integer, dimension(2**get_num_dims(self))  :: nodes
                  integer                                    :: nsize

                  nsize = 2**get_num_dims(self)
                  call get_mesh_cell_corner_nodes(self%this,cell,nodes,nsize)

              end function get_cell_corner_nodes

! Return an array of the nodes that comprise a cell face.
              function get_cell_face_nodes(self, cell, face) result(nodes)
                  type(CAR_CU_Mesh), intent(in)                   :: self
                  integer, intent(in)                             :: cell, face
                  integer, dimension(2**(get_num_dims(self)-1))   :: nodes
                  integer                                         :: nsize

                  nsize = 2**(get_num_dims(self) - 1)
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
                                  get_num_dims(self))      :: vertices
                  real*8, dimension(get_num_nodes(self) *               &
                                  get_num_dims(self))      :: ret_vert
                  integer                                  :: node, dir, nsize

                  nsize = get_num_nodes(self) * get_num_dims(self)
                  call get_mesh_vertices(self%this, ret_vert, nsize)

                  node = 1
                  do while (node .le. get_num_nodes(self))
                      dir = 1
                      do while (dir .le. get_num_dims(self))
                          vertices(node, dir) =                         &
                              ret_vert(get_num_dims(self)*(node-1) + dir)
                          dir = dir + 1
                      end do
                      node = node + 1
                  end do

              end function get_vertices

! Return an array containing the vertices for all of the cell corner nodes.

              function get_corner_node_vertices(self) result(vertices)
                  type(CAR_CU_Mesh), intent(in)             :: self
                  real*8, dimension(get_num_corner_nodes(self),         &
                                    get_num_dims(self))     :: vertices
                  real*8, dimension(get_num_corner_nodes(self) *        &
                                    get_num_dims(self))     :: ret_vert
                  integer                                   :: node, dir, nsize

                  nsize = get_num_corner_nodes(self) * get_num_dims(self)
                  call get_mesh_corner_node_vertices(self%this,ret_vert,nsize)

                  node = 1
                  do while (node .le. get_num_corner_nodes(self))
                      dir = 1
                      do while (dir .le. get_num_dims(self))
                          vertices(node, dir) =                         &
                              ret_vert(get_num_dims(self)*(node-1) + dir)
                          dir = dir + 1
                      end do
                      node = node + 1
                  end do

              end function get_corner_node_vertices

! Return an array containing the vertices for all of the face_centered nodes.

              function get_face_centered_node_vertices(self) result(vertices)
                  type(CAR_CU_Mesh), intent(in)            :: self
                  real*8, dimension(get_num_face_nodes(self),           &
                                    get_num_dims(self))    :: vertices
                  real*8, dimension(get_num_face_nodes(self) *          &
                                    get_num_dims(self))    :: ret_vert
                  integer                                  :: node, dir, nsize

                  nsize = get_num_face_nodes(self) * get_num_dims(self)
                  call get_mesh_face_cen_node_vertices(self%this,       &
                                                       ret_vert, nsize)

                  node = 1
                  do while (node .le. get_num_face_nodes(self))
                      dir = 1
                      do while (dir .le. get_num_dims(self))
                          vertices(node, dir) =                         &
                              ret_vert(get_num_dims(self)*(node-1) + dir)
                          dir = dir + 1
                      end do
                      node = node + 1
                  end do

              end function get_face_centered_node_vertices

! Return an array with all of a cell's vertices
              function get_cell_vertices(self, cell) result(vertices)
                  type(CAR_CU_Mesh), intent(in)             :: self
                  integer, intent(in)                       :: cell
                  real*8, dimension(2**get_num_dims(self),              &
                                    get_num_dims(self))     :: vertices
                  real*8, dimension(2**get_num_dims(self) *             &
                                    get_num_dims(self))     :: ret_vert

                  integer                                   :: node, nsize, dir

                  nsize = 2**get_num_dims(self) * get_num_dims(self)
                  call get_mesh_cell_vertices(self%this, cell,ret_vert, nsize)

                  node = 1
                  do while (node .le. 2**get_num_dims(self))
                      dir = 1
                      do while (dir .le. get_num_dims(self))
                          vertices(node, dir) =                         &
                              ret_vert(get_num_dims(self)*(node-1) + dir)
                          dir = dir + 1
                      end do
                      node = node + 1
                  end do

              end function get_cell_vertices

! Return an array of all of the cell face vertices
              function get_cell_face_vertices(self,cell,face) result(vertices)
                  type(CAR_CU_Mesh), intent(in)             :: self
                  integer, intent(in)                       :: cell, face
                  real*8, dimension(2*(get_num_dims(self)-1),           &
                                    get_num_dims(self))     :: vertices
                  real*8, dimension(2*(get_num_dims(self)-1) *          &
                                    get_num_dims(self))     :: ret_vert
                  integer                                   :: node, nsize, dir

                  nsize = 2 * (get_num_dims(self) - 1) * get_num_dims(self)
                  call get_mesh_cell_face_vertices(self%this, cell,     &
                                                   face, ret_vert, nsize)
                  node = 1
                  do while (node .le. 2*(get_num_dims(self) - 1))
                      dir = 1
                      do while (dir .le. get_num_dims(self))
                          vertices(node, dir) =                         &
                              ret_vert(get_num_dims(self)*(node-1) + dir)
                          dir = dir + 1
                      end do
                      node = node + 1
                  end do

              end function get_cell_face_vertices

! Return a single node's vertices
              function get_node_vertices(self, node) result(vertices)
                  type(CAR_CU_Mesh), intent(in)            :: self
                  integer, intent(in)                      :: node
                  real*8, dimension(get_num_dims(self))    :: vertices
                  integer                                  :: nsize

                  nsize = get_num_dims(self)
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
              function get_integer_CCSF_all(self)           result(data)
                  type(integer_CCSF), intent(in)                :: self
                  integer, dimension(get_num_cells(self%mesh))  :: data
                  integer                                       :: data_size

                  data_size = get_num_cells(self%mesh)
                  call get_mesh_ccsf_i(self%mesh%this, self%this, data, &
                                       data_size)

              end function get_integer_CCSF_all

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
              subroutine set_integer_CCSF_all(self, data)
                  type(integer_CCSF), intent(in)                :: self
                  integer, intent(in),                                  &
                           dimension(get_num_cells(self%mesh))  :: data
                  integer                                       :: data_size

                  data_size = get_num_cells(self%mesh)
                  call set_mesh_ccsf_i(self%mesh%this, self%this, data, &
                                       data_size)

              end subroutine set_integer_CCSF_all

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
              function get_real_CCSF_all(self)              result(data)
                  type(real_CCSF), intent(in)                   :: self
                  real*8, dimension(get_num_cells(self%mesh))   :: data
                  integer                                       :: data_size

                  data_size = get_num_cells(self%mesh)
                  call get_mesh_ccsf_d(self%mesh%this, self%this, data, &
                                       data_size)

              end function get_real_CCSF_all

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
              subroutine set_real_CCSF_all(self, data)
                  type(real_CCSF), intent(in)                   :: self
                  real*8, intent(in),                                   &
                          dimension(get_num_cells(self%mesh))   :: data
                  integer                                       :: data_size

                  data_size = get_num_cells(self%mesh)
                  call set_mesh_ccsf_d(self%mesh%this, self%this, data, &
                                       data_size)

              end subroutine set_real_CCSF_all

! Set a cell value from a C++ CAR_CU_Mesh double CCSF class object (self).
              subroutine set_real_CCSF_cell(self, cell, data)
                  type(real_CCSF), intent(in)              :: self
                  integer, intent(in)                      :: cell
                  real*8, intent(in)                       :: data

                  call set_mesh_ccsf_d_cell(self%mesh%this, self%this,  &
                                            cell, data)

              end subroutine set_real_CCSF_cell

!===========================================================================
! integer CCVF class objects
!===========================================================================
! Return an entire C++ CAR_CU_Mesh integer CCVF class object (self).
              function get_integer_CCVF_all(self)             result(data)
                  type(integer_CCVF), intent(in)               :: self
                  integer, dimension(get_num_cells(self%mesh),          &
                                     self%vec_size)            :: data
                  integer, dimension(get_num_cells(self%mesh) *         &
                                     self%vec_size)            :: ret_data
                  integer                                      :: data_size
                  integer                                      :: cell,dim

                  data_size = get_num_cells(self%mesh) * self%vec_size
                  call get_mesh_ccvf_i(self%mesh%this, self%this,       &
                                       ret_data, data_size)

                  cell = 1
                  do while(cell .le. get_num_cells(self%mesh))
                      dim = 1
                      do while (dim .le. self%vec_size)
                          data(cell, dim) = ret_data(self%vec_size *    &
                              (cell-1) + dim)
                          dim = dim + 1
                      end do
                      cell = cell + 1
                  end do

              end function get_integer_CCVF_all

! Return all of the dim values for a cell in a C++ CAR_CU_Mesh integer CCVF 
! class object (self).
              function get_integer_CCVF_cell(self, cell)      result(data)
                  type(integer_CCVF), intent(in)               :: self
                  integer, intent(in)                          :: cell
                  integer, dimension(self%vec_size)            :: data
                  integer                                      :: data_size

                  data_size = self%vec_size
                  call get_mesh_ccvf_i_cell(self%mesh%this, self%this,  &
                                            cell, data, data_size)

              end function get_integer_CCVF_cell

! Return a cell dim value from a C++ CAR_CU_Mesh integer CCVF class object 
! (self).
              function get_integer_CCVF_cell_dim(self, cell, dim)       &
                                                     result(data)
                  type(integer_CCVF), intent(in)         :: self
                  integer, intent(in)                    :: cell, dim
                  integer                                :: data

                  call get_mesh_ccvf_i_cell_dim(self%mesh%this,         &
                                    self%this, cell, dim, data)

              end function get_integer_CCVF_cell_dim

! Set an entire C++ CAR_CU_Mesh integer CCVF class object (self) (can also
! be done at initialization using the constructor).
              subroutine set_integer_CCVF_all(self, data)
                  type(integer_CCVF), intent(in)               :: self
                  integer, intent(in),                                  &
                           dimension(get_num_cells(self%mesh),          &
                                     self%vec_size)            :: data
                  integer, dimension(get_num_cells(self%mesh) *         &
                                     self%vec_size)            :: ret_data
                  integer                                      :: data_size
                  integer                                      :: cell,dim

                  cell = 1
                  do while(cell .le. get_num_cells(self%mesh))
                      dim = 1
                      do while (dim .le. self%vec_size)
                          ret_data(self%vec_size * (cell - 1) + dim) =  &
                              data(cell, dim)
                          dim = dim + 1
                      end do
                      cell = cell + 1
                  end do

                  data_size = get_num_cells(self%mesh) * self%vec_size
                  call set_mesh_ccvf_i(self%mesh%this, self%this,       &
                                        ret_data, data_size)

              end subroutine set_integer_CCVF_all

! Set all of the dim values for a cell in a C++ CAR_CU_Mesh integer CCVF 
! class object (self).
              subroutine set_integer_CCVF_cell(self, cell, data)
                  type(integer_CCVF), intent(in)               :: self
                  integer, intent(in)                          :: cell
                  integer, intent(in),                                  &
                           dimension(self%vec_size)            :: data
                  integer                                      :: data_size

                  data_size = self%vec_size
                  call set_mesh_ccvf_i_cell(self%mesh%this, self%this,  &
                                             cell, data, data_size)

              end subroutine set_integer_CCVF_cell

! Set a cell dim value for a C++ CAR_CU_Mesh integer CCVF class object 
! (self).
              subroutine set_integer_CCVF_cell_dim(self, cell, dim, data)
                  type(integer_CCVF), intent(in)           :: self
                  integer, intent(in)                      :: cell, dim
                  integer, intent(in)                      :: data

                  call set_mesh_ccvf_i_cell_dim(self%mesh%this,         &
                                     self%this, cell, dim, data)

              end subroutine set_integer_CCVF_cell_dim

!===========================================================================
! double CCVF class objects
!===========================================================================
! Return an entire C++ CAR_CU_Mesh double CCVF class object (self).
              function get_real_CCVF_all(self)            result(data)
                  type(real_CCVF), intent(in)                 :: self
                  real*8, dimension(get_num_cells(self%mesh),           &
                                    self%vec_size)            :: data
                  real*8, dimension(get_num_cells(self%mesh) *          &
                                    self%vec_size)            :: ret_data
                  integer                                     :: data_size
                  integer                                     :: cell,dim

                  data_size = get_num_cells(self%mesh) *  self%vec_size
                  call get_mesh_ccvf_d(self%mesh%this, self%this,       &
                                       ret_data, data_size)

                  cell = 1
                  do while(cell .le. get_num_cells(self%mesh))
                      dim = 1
                      do while (dim .le. self%vec_size)
                          data(cell, dim) =                             &
                              ret_data(self%vec_size * (cell-1) + dim)
                          dim = dim + 1
                      end do
                      cell = cell + 1
                  end do

              end function get_real_CCVF_all

! Return all of the dim values for a cell in a C++ CAR_CU_Mesh double CCVF 
! class object (self).
              function get_real_CCVF_cell(self, cell)     result(data)
                  type(real_CCVF), intent(in)                :: self
                  integer, intent(in)                        :: cell
                  real*8, dimension(self%vec_size)           :: data
                  integer                                    :: data_size

                  data_size = self%vec_size
                  call get_mesh_ccvf_d_cell(self%mesh%this, self%this,  &
                                            cell, data, data_size)

              end function get_real_CCVF_cell

! Return a cell dim value from a C++ CAR_CU_Mesh double CCVF class object
! (self).
              function get_real_CCVF_cell_dim(self, cell, dim) result(data)
                  type(real_CCVF), intent(in)            :: self
                  integer, intent(in)                    :: cell, dim
                  real*8                                 :: data

                  call get_mesh_ccvf_d_cell_dim(self%mesh%this,         &
                                     self%this, cell, dim, data)

              end function get_real_CCVF_cell_dim

! Set an entire C++ CAR_CU_Mesh double CCVF class object (self) (can also
! be done at initialization using the constructor).
              subroutine set_real_CCVF_all(self, data)
                  type(real_CCVF), intent(in)                 :: self
                  real*8, intent(in),                                   &
                          dimension(get_num_cells(self%mesh),           &
                                    self%vec_size)            :: data
                  real*8, dimension(get_num_cells(self%mesh) *          &
                                    self%vec_size)            :: ret_data
                  integer                                     :: data_size
                  integer                                     :: cell,dim

                  cell = 1
                  do while(cell .le. get_num_cells(self%mesh))
                      dim = 1
                      do while (dim .le. self%vec_size)
                          ret_data(self%vec_size * (cell-1) + dim) =    &
                              data(cell, dim)
                          dim = dim + 1
                      end do
                      cell = cell + 1
                  end do

                  data_size = get_num_cells(self%mesh) * self%vec_size
                  call set_mesh_ccvf_d(self%mesh%this, self%this,       &
                                       ret_data, data_size)

              end subroutine set_real_CCVF_all

! Set all of the dim values for a cell in a C++ CAR_CU_Mesh double CCVF class
! object (self).
              subroutine set_real_CCVF_cell(self, cell, data)
                  type(real_CCVF), intent(in)                 :: self
                  integer, intent(in)                         :: cell
                  real*8, intent(in),                                   &
                          dimension(self%vec_size)            :: data
                  integer                                     :: data_size

                  data_size = self%vec_size
                  call set_mesh_ccvf_d_cell(self%mesh%this, self%this,  &
                                            cell, data, data_size)

              end subroutine set_real_CCVF_cell

! Set a cell dim value from a C++ CAR_CU_Mesh double CCVF class object 
! (self).
              subroutine set_real_CCVF_cell_dim(self, cell, dim, data)
                  type(real_CCVF), intent(in)              :: self
                  integer, intent(in)                      :: cell, dim
                  real*8, intent(in)                       :: data

                  call set_mesh_ccvf_d_cell_dim(self%mesh%this,        &
                                    self%this, cell, dim, data)

              end subroutine set_real_CCVF_cell_dim

!===========================================================================
! integer FCSF class objects
!===========================================================================
! Return an entire C++ CAR_CU_Mesh integer FCSF class object (self).
              function get_integer_FCSF_all(self)           result(data)
                  type(integer_FCSF), intent(in)                :: self
                  integer, dimension(get_num_face_nodes(self%mesh)) :: data
                  integer                                       :: data_size

                  data_size = get_num_face_nodes(self%mesh)
                  call get_mesh_fcsf_i(self%mesh%this, self%this, data, &
                                       data_size)

              end function get_integer_FCSF_all

! Return all of the face values for a cell in a C++ CAR_CU_Mesh integer FCSF 
! class object (self).
              function get_integer_FCSF_cell(self, cell)          result(data)
                  type(integer_FCSF), intent(in)                   :: self
                  integer, intent(in)                              :: cell
                  integer, dimension(2 * get_num_dims(self%mesh))  :: data
                  integer                                          :: data_size

                  data_size = 2 * get_num_dims(self%mesh)
                  call get_mesh_fcsf_i_cell(self%mesh%this, self%this,  &
                                            cell, data, data_size)

              end function get_integer_FCSF_cell

! Return a cell face value from a C++ CAR_CU_Mesh integer FCSF class object 
! (self).
              function get_integer_FCSF_cell_face(self, cell, face)     &
                                                       result(data)
                  type(integer_FCSF), intent(in)           :: self
                  integer, intent(in)                      :: cell, face
                  integer                                  :: data

                  call get_mesh_fcsf_i_cell_face(self%mesh%this,        &
                                    self%this, cell, face, data)

              end function get_integer_FCSF_cell_face

! Set an entire C++ CAR_CU_Mesh integer FCSF class object (self) (can also
! be done at initialization using the constructor).
              subroutine set_integer_FCSF_all(self, data)
                  type(integer_FCSF), intent(in)                   :: self
                  integer, intent(in),                                  &
                           dimension(get_num_face_nodes(self%mesh)):: data
                  integer                                          :: data_size

                  data_size = get_num_face_nodes(self%mesh)
                  call set_mesh_fcsf_i(self%mesh%this, self%this, data, &
                                       data_size)

              end subroutine set_integer_FCSF_all

! Set all of the face values for a cell in a C++ CAR_CU_Mesh integer FCSF class
! object (self).
              subroutine set_integer_FCSF_cell(self, cell, data)
                  type(integer_FCSF), intent(in)                   :: self
                  integer, intent(in)                              :: cell
                  integer, intent(in),                                  &
                           dimension(2 * get_num_dims(self%mesh))  :: data
                  integer                                          :: data_size

                  data_size = 2 * get_num_dims(self%mesh)
                  call set_mesh_fcsf_i_cell(self%mesh%this, self%this,  &
                                            cell, data, data_size)

              end subroutine set_integer_FCSF_cell

! Set a cell face value for a C++ CAR_CU_Mesh integer FCSF class object (self).
              subroutine set_integer_FCSF_cell_face(self, cell, face, data)
                  type(integer_FCSF), intent(in)           :: self
                  integer, intent(in)                      :: cell, face
                  integer, intent(in)                      :: data

                  call set_mesh_fcsf_i_cell_face(self%mesh%this,        &
                                    self%this, cell, face, data)

              end subroutine set_integer_FCSF_cell_face

!===========================================================================
! double FCSF class objects
!===========================================================================
! Return an entire C++ CAR_CU_Mesh double FCSF class object (self).
              function get_real_FCSF_all(self)                 result(data)
                  type(real_FCSF), intent(in)                      :: self
                  real*8, dimension(get_num_face_nodes(self%mesh)) :: data
                  integer                                          :: data_size

                  data_size = get_num_face_nodes(self%mesh)
                  call get_mesh_fcsf_d(self%mesh%this, self%this, data, &
                                       data_size)

              end function get_real_FCSF_all

! Return all of the face values for a cell in a C++ CAR_CU_Mesh double FCSF 
! class object (self).
              function get_real_FCSF_cell(self, cell)          result(data)
                  type(real_FCSF), intent(in)                      :: self
                  integer, intent(in)                              :: cell
                  real*8, dimension(2 * get_num_dims(self%mesh))   :: data
                  integer                                          :: data_size

                  data_size = 2 * get_num_dims(self%mesh)
                  call get_mesh_fcsf_d_cell(self%mesh%this, self%this,  &
                                            cell, data, data_size)

              end function get_real_FCSF_cell

! Return a cell face value from a C++ CAR_CU_Mesh double FCSF class object
! (self).
              function get_real_FCSF_cell_face(self, cell, face) result(data)
                  type(real_FCSF), intent(in)             :: self
                  integer, intent(in)                     :: cell, face
                  real*8                                  :: data

                  call get_mesh_fcsf_d_cell_face(self%mesh%this,        &
                                    self%this, cell, face, data)

              end function get_real_FCSF_cell_face

! Set an entire C++ CAR_CU_Mesh double FCSF class object (self) (can also
! be done at initialization using the constructor).
              subroutine set_real_FCSF_all(self, data)
                  type(real_FCSF), intent(in)                   :: self
                  real*8, intent(in),                                   &
                          dimension(get_num_face_nodes(self%mesh)) :: data
                  integer                                       :: data_size

                  data_size = get_num_face_nodes(self%mesh)
                  call set_mesh_fcsf_d(self%mesh%this, self%this, data, &
                                       data_size)

              end subroutine set_real_FCSF_all

! Set all of the face values for a cell in a C++ CAR_CU_Mesh double FCSF class
! object (self).
              subroutine set_real_FCSF_cell(self, cell, data)
                  type(real_FCSF), intent(in)                      :: self
                  integer, intent(in)                              :: cell
                  real*8, intent(in),                                   &
                          dimension(2 * get_num_dims(self%mesh))   :: data
                  integer                                          :: data_size

                  data_size = 2 * get_num_dims(self%mesh)
                  call set_mesh_fcsf_d_cell(self%mesh%this, self%this,  &
                                            cell, data, data_size)

              end subroutine set_real_FCSF_cell

! Set a cell face value from a C++ CAR_CU_Mesh double FCSF class object (self).
              subroutine set_real_FCSF_cell_face(self, cell, face, data)
                  type(real_FCSF), intent(in)              :: self
                  integer, intent(in)                      :: cell, face
                  real*8, intent(in)                       :: data

                  call set_mesh_fcsf_d_cell_face(self%mesh%this,        &
                                    self%this, cell, face, data)

              end subroutine set_real_FCSF_cell_face

!===========================================================================
! integer FCDSF class objects
!===========================================================================
! Return an entire C++ CAR_CU_Mesh integer FCDSF class object (self).
              function get_integer_FCDSF_all(self)             result(data)
                  type(integer_FCDSF), intent(in)                  :: self
                  integer, dimension(get_num_cells(self%mesh),          &
                                     2 * get_num_dims(self%mesh))  :: data
                  integer, dimension(get_num_cells(self%mesh) *         &
                                     2 * get_num_dims(self%mesh))  :: ret_data
                  integer                                          :: data_size
                  integer                                          :: cell,face

                  data_size = get_num_cells(self%mesh) *                &
                              2 * get_num_dims(self%mesh)
                  call get_mesh_fcdsf_i(self%mesh%this, self%this,      &
                                        ret_data, data_size)

                  cell = 1
                  do while(cell .le. get_num_cells(self%mesh))
                      face = 1
                      do while (face .le. 2 * get_num_dims(self%mesh))
                          data(cell, face) =                            &
                              ret_data(2 * get_num_dims(self%mesh) *    &
                              (cell-1) + face)
                          face = face + 1
                      end do
                      cell = cell + 1
                  end do

              end function get_integer_FCDSF_all

! Return all of the face values for a cell in a C++ CAR_CU_Mesh integer FCDSF 
! class object (self).
              function get_integer_FCDSF_cell(self, cell)      result(data)
                  type(integer_FCDSF), intent(in)                  :: self
                  integer, intent(in)                              :: cell
                  integer, dimension(2 * get_num_dims(self%mesh))  :: data
                  integer                                          :: data_size

                  data_size = 2 * get_num_dims(self%mesh)
                  call get_mesh_fcdsf_i_cell(self%mesh%this, self%this, &
                                             cell, data, data_size)

              end function get_integer_FCDSF_cell

! Return a cell face value from a C++ CAR_CU_Mesh integer FCDSF class object 
! (self).
              function get_integer_FCDSF_cell_face(self, cell, face)    &
                                                       result(data)
                  type(integer_FCDSF), intent(in)          :: self
                  integer, intent(in)                      :: cell, face
                  integer                                  :: data

                  call get_mesh_fcdsf_i_cell_face(self%mesh%this,       &
                                     self%this, cell, face, data)

              end function get_integer_FCDSF_cell_face

! Set an entire C++ CAR_CU_Mesh integer FCDSF class object (self) (can also
! be done at initialization using the constructor).
              subroutine set_integer_FCDSF_all(self, data)
                  type(integer_FCDSF), intent(in)                  :: self
                  integer, intent(in),                                  &
                           dimension(get_num_cells(self%mesh),          &
                                     2 * get_num_dims(self%mesh))  :: data
                  integer, dimension(get_num_cells(self%mesh) *         &
                                     2 * get_num_dims(self%mesh))  :: ret_data
                  integer                                          :: data_size
                  integer                                          :: cell,face

                  cell = 1
                  do while(cell .le. get_num_cells(self%mesh))
                      face = 1
                      do while (face .le. 2 * get_num_dims(self%mesh))
                          ret_data(2 * get_num_dims(self%mesh) *        &
                              (cell - 1) + face) = data(cell, face)
                          face = face + 1
                      end do
                      cell = cell + 1
                  end do

                  data_size = get_num_cells(self%mesh) *                &
                              2 * get_num_dims(self%mesh)
                  call set_mesh_fcdsf_i(self%mesh%this, self%this,      &
                                        ret_data, data_size)

              end subroutine set_integer_FCDSF_all

! Set all of the face values for a cell in a C++ CAR_CU_Mesh integer FCDSF 
! class object (self).
              subroutine set_integer_FCDSF_cell(self, cell, data)
                  type(integer_FCDSF), intent(in)                  :: self
                  integer, intent(in)                              :: cell
                  integer, intent(in),                                  &
                           dimension(2 * get_num_dims(self%mesh))  :: data
                  integer                                          :: data_size

                  data_size = 2 * get_num_dims(self%mesh)
                  call set_mesh_fcdsf_i_cell(self%mesh%this, self%this, &
                                             cell, data, data_size)

              end subroutine set_integer_FCDSF_cell

! Set a cell face value for a C++ CAR_CU_Mesh integer FCDSF class object 
! (self).
              subroutine set_integer_FCDSF_cell_face(self, cell, face, data)
                  type(integer_FCDSF), intent(in)          :: self
                  integer, intent(in)                      :: cell, face
                  integer, intent(in)                      :: data

                  call set_mesh_fcdsf_i_cell_face(self%mesh%this,       &
                                     self%this, cell, face, data)

              end subroutine set_integer_FCDSF_cell_face

!===========================================================================
! double FCDSF class objects
!===========================================================================
! Return an entire C++ CAR_CU_Mesh double FCDSF class object (self).
              function get_real_FCDSF_all(self)               result(data)
                  type(real_FCDSF), intent(in)                    :: self
                  real*8, dimension(get_num_cells(self%mesh),           &
                                    2 * get_num_dims(self%mesh))  :: data
                  real*8, dimension(get_num_cells(self%mesh) *          &
                                    2 * get_num_dims(self%mesh))  :: ret_data
                  integer                                         :: data_size
                  integer                                         :: cell,face

                  data_size = get_num_cells(self%mesh) *                &
                              2 * get_num_dims(self%mesh)
                  call get_mesh_fcdsf_d(self%mesh%this, self%this,      &
                                        ret_data, data_size)

                  cell = 1
                  do while(cell .le. get_num_cells(self%mesh))
                      face = 1
                      do while (face .le. 2 * get_num_dims(self%mesh))
                          data(cell, face) =                            &
                              ret_data(2 * get_num_dims(self%mesh) *    &
                              (cell-1) + face)
                          face = face + 1
                      end do
                      cell = cell + 1
                  end do

              end function get_real_FCDSF_all

! Return all of the face values for a cell in a C++ CAR_CU_Mesh double FCDSF 
! class object (self).
              function get_real_FCDSF_cell(self, cell)        result(data)
                  type(real_FCDSF), intent(in)                    :: self
                  integer, intent(in)                             :: cell
                  real*8, dimension(2 * get_num_dims(self%mesh))  :: data
                  integer                                         :: data_size

                  data_size = 2 * get_num_dims(self%mesh)
                  call get_mesh_fcdsf_d_cell(self%mesh%this, self%this, &
                                             cell, data, data_size)

              end function get_real_FCDSF_cell

! Return a cell face value from a C++ CAR_CU_Mesh double FCDSF class object
! (self).
              function get_real_FCDSF_cell_face(self, cell, face) result(data)
                  type(real_FCDSF), intent(in)            :: self
                  integer, intent(in)                     :: cell, face
                  real*8                                  :: data

                  call get_mesh_fcdsf_d_cell_face(self%mesh%this,       &
                                     self%this, cell, face, data)

              end function get_real_FCDSF_cell_face

! Set an entire C++ CAR_CU_Mesh double FCDSF class object (self) (can also
! be done at initialization using the constructor).
              subroutine set_real_FCDSF_all(self, data)
                  type(real_FCDSF), intent(in)                    :: self
                  real*8, intent(in),                                   &
                          dimension(get_num_cells(self%mesh),           &
                                    2 * get_num_dims(self%mesh))  :: data
                  real*8, dimension(get_num_cells(self%mesh) *          &
                                    2 * get_num_dims(self%mesh))  :: ret_data
                  integer                                         :: data_size
                  integer                                         :: cell,face

                  cell = 1
                  do while(cell .le. get_num_cells(self%mesh))
                      face = 1
                      do while (face .le. 2 * get_num_dims(self%mesh))
                          ret_data(2 * get_num_dims(self%mesh) *        &
                              (cell-1) + face) = data(cell, face)
                          face = face + 1
                      end do
                      cell = cell + 1
                  end do

                  data_size = get_num_cells(self%mesh) *                &
                              2 * get_num_dims(self%mesh)
                  call set_mesh_fcdsf_d(self%mesh%this, self%this,      &
                                        ret_data, data_size)

              end subroutine set_real_FCDSF_all

! Set all of the face values for a cell in a C++ CAR_CU_Mesh double FCDSF class
! object (self).
              subroutine set_real_FCDSF_cell(self, cell, data)
                  type(real_FCDSF), intent(in)                     :: self
                  integer, intent(in)                              :: cell
                  real*8, intent(in),                                   &
                          dimension(2 * get_num_dims(self%mesh))   :: data
                  integer                                          :: data_size

                  data_size = 2 * get_num_dims(self%mesh)
                  call set_mesh_fcdsf_d_cell(self%mesh%this, self%this, &
                                             cell, data, data_size)

              end subroutine set_real_FCDSF_cell

! Set a cell face value from a C++ CAR_CU_Mesh double FCDSF class object 
! (self).
              subroutine set_real_FCDSF_cell_face(self, cell, face, data)
                  type(real_FCDSF), intent(in)              :: self
                  integer, intent(in)                       :: cell, face
                  real*8, intent(in)                        :: data

                  call set_mesh_fcdsf_d_cell_face(self%mesh%this,       &
                                     self%this, cell, face, data)

              end subroutine set_real_FCDSF_cell_face

      end module CAR_CU_Mesh_Class


!---------------------------------------------------------------------------
!                              end of amr_mesh/Shadow_CAR_CU_Mesh.f90
!---------------------------------------------------------------------------
