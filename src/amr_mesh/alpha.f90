!----------------------------------*-F90-*----------------------------------
! alpha.f90
! B.T. Adams (bta@lanl.gov)
! 30 Sept 99
!---------------------------------------------------------------------------
! @> alpha program file
!---------------------------------------------------------------------------
!===========================================================================
! alpha - 
!
! Purpose : Main driver program for the zathros/centauri combination.
!
! revision history:
! -----------------
!  0) original
! 
!===========================================================================

      program alpha
!===========================================================================
! Shadow Interface modules - have to have this stuff
!===========================================================================

          USE CAR_CU_Interface_Class
          USE CAR_CU_RTT_Format_Class
          USE CAR_CU_Builder_Class
          USE CAR_CU_Mesh_Class

          implicit none
          integer narg, iargc, fnlgth
          type(CAR_CU_Interface)  :: interface_class
          type(CAR_CU_RTT_Format) :: rtt_format_class
          type(CAR_CU_Builder)    :: builder_class
          type(CAR_CU_Mesh)       :: mesh_class

!===========================================================================
! Define variables needed just for testing the shadow interface functions
!===========================================================================

          type(integer_CCSF)      ::   int_CCSF_1,   int_CCSF_2
          type(real_CCSF)         ::  real_CCSF_1,  real_CCSF_2

          type(integer_CCVF)      ::   int_CCVF_1,   int_CCVF_2,        &
                                       int_CCVF_3,   int_CCVF_4
          type(real_CCVF)         ::  real_CCVF_1,  real_CCVF_2,        &
                                      real_CCVF_3,  real_CCVF_4

          type(integer_FCSF)      ::   int_FCSF_1,   int_FCSF_2
          type(real_FCSF)         ::  real_FCSF_1,  real_FCSF_2

          type(integer_FCDSF)     ::  int_FCDSF_1,  int_FCDSF_2
          type(real_FCDSF)        :: real_FCDSF_1, real_FCDSF_2

          type(integer_FCVF)      ::   int_FCVF_1,   int_FCVF_2,        &
                                       int_FCVF_3,   int_FCVF_4
          type(real_FCVF)         ::  real_FCVF_1,  real_FCVF_2,        &
                                      real_FCVF_3,  real_FCVF_4

          type(integer_NCSF)      ::   int_NCSF_1,   int_NCSF_2,        &
                                       int_NCSF_3,   int_NCSF_4
          type(real_NCSF)         ::  real_NCSF_1,  real_NCSF_2,        &
                                      real_NCSF_3,  real_NCSF_4
          type(integer_NCVF)      ::   int_NCVF_1,   int_NCVF_2,        &
                                       int_NCVF_3,   int_NCVF_4,        &
                                       int_NCVF_5,   int_NCVF_6
          type(real_NCVF)         ::  real_NCVF_1,  real_NCVF_2,        &
                                      real_NCVF_3,  real_NCVF_4,        &
                                      real_NCVF_5,  real_NCVF_6

          integer ndims, dir, ncells, nnodes, ncnodes, nfnodes, cell,   &
              face, node

          integer, dimension (:,:), allocatable :: num_adj_cell_faces,  &
              adj_cells_faces, cells_nodes, dis_face_generation,        &
              cells_dirs, corner_nodes_faces_gen, nodes_dirs,           &
              corner_nodes_dirs, face_nodes_dirs, face_nodes_faces_gen

          real*8, dimension (:,:), allocatable  :: nodes_vertices,      &
              corner_nodes_vertices, centered_nodes_vertices,           &
              cells_faces_area, face_nodes_vertices, cells_min_coords,  &
              cells_mid_coords, cells_max_coords, cells_dirs_width,     &
              cell_nodes_coords, dis_face_nodes_area,                   &
              corner_nodes_faces_area, face_nodes_faces_area

          integer, dimension (:), allocatable   :: generation,          &
              cell_face_nodes, cell_faces_centered_node,                &
              cell_faces_specific_node, cell_nodes, cell_corner_nodes,  &
              cell_face_centered_nodes, face_generation,                &
              nodes_generation, corner_nodes_generation, dirs

          real*8, dimension (:), allocatable    :: volume,              &
              cell_faces_area, mesh_min_coords, mesh_max_coords,        &
              cell_nodes_vertices, face_nodes_area, nodes_area,         &
              corner_nodes_area

!===========================================================================
! Input the command line arguments - input file name followed by anything to
! activate the verbose switch (typically v or the word verbose).
!===========================================================================

          narg = iargc()
          if (narg .gt. 0) then
              call getarg(1,interface_class%file_name)
              ! Have to add a terminating null character to the end of the 
              ! file name for conversion to a C++ string object
              fnlgth = len_trim(interface_class%file_name)
              interface_class%file_name(fnlgth + 1:fnlgth + 1) =  ACHAR(0)
              if (narg .gt. 1) then
                  interface_class%verbose = .TRUE.
              else
                  interface_class%verbose = .FALSE.
              end if
          else
              stop
          end if

!===========================================================================
! Create a C++ CAR_CU_Interface class object. This also constructs a C++ 
! RTT_Format class object and parses both the input deck and the RTT Format 
! mesh file specified therein. The addresses of both the new CAR_CU_Interace 
! and RTT_Format class objects are set.
!===========================================================================

          call construct_Interface(interface_class)
          rtt_format_class%this = interface_class%rtt_format

!===========================================================================
! Create a C++ CAR_CU_Builder class object. This also constructs the C++ 
! Coord_sys, Layout, and CAR_CU_Mesh class objects. The addresses of both 
! the new CAR_CU_Builder and CAR_CU_Mesh class objects are set.
!===========================================================================

          call construct_Builder(builder_class, interface_class)
          mesh_class%this = builder_class%mesh

!===========================================================================
! Get rid of RTT_Format, CAR_CU_Interface, and CAR_CU_Builder class objects 
! since they are no longer needed. 
!===========================================================================

          call destruct_RTT_Format(rtt_format_class)
          call destruct_Interface(interface_class)
          call destruct_Builder(builder_class)

!===========================================================================
! Test shadowed mesh accessor functions
!===========================================================================

          ! Test scalar values
          ndims = get_num_dims(mesh_class)
          ncells = get_num_cells(mesh_class)
          nnodes = get_num_nodes(mesh_class)
          ncnodes = get_num_corner_nodes(mesh_class)
          nfnodes = get_num_face_nodes(mesh_class)

          ! Allocate memory for test variable arrays
          ! Cell-centered values
          allocate(volume(ncells))
          allocate(generation(ncells))
          ! Face-dependent cell values
          allocate(num_adj_cell_faces(ncells, 2 * ndims))
          allocate(adj_cells_faces(ncells, 2 * ndims))
          allocate(cells_faces_area(ncells, 2 * ndims))
          allocate(cell_faces_centered_node( 2 * ndims))
          ! Face-dependent values
          allocate(cell_faces_area( 2 * ndims))
          ! Direction-dependent values
          allocate(mesh_min_coords(ndims))
          allocate(mesh_max_coords(ndims))
          allocate(dirs(ndims))
          ! Direction-dependent cell values
          allocate(cells_min_coords(ncells, ndims))
          allocate(cells_mid_coords(ncells, ndims))
          allocate(cells_max_coords(ncells, ndims))
          allocate(cells_dirs_width(ncells, ndims))
          allocate(cells_dirs(ncells, ndims))
          ! Node-dependent cell values
          allocate(cell_face_nodes(2 * (ndims - 1)))
          allocate(cell_faces_specific_node(2 * ndims))
          allocate(cell_nodes(2**ndims + 2*ndims))
          allocate(cell_corner_nodes(2**ndims))
          allocate(cell_face_centered_nodes(2 * ndims))
          allocate(cells_nodes(ncells, (2**ndims + 2*ndims)))
          ! Node-centered values
          allocate(nodes_vertices(nnodes, ndims))
          allocate(corner_nodes_vertices(ncnodes, ndims))
          allocate(centered_nodes_vertices(nfnodes, ndims))
          allocate(face_nodes_vertices((2*(ndims - 1)), ndims))
          allocate(cell_nodes_coords((2**ndims), ndims))
          allocate(cell_nodes_vertices(2**ndims + 2*ndims))
          allocate(face_generation(nfnodes))
          allocate(face_nodes_area(nfnodes))
          allocate(dis_face_generation(ncells, 2 * ndims))
          allocate(dis_face_nodes_area(ncells, 2 * ndims))
          allocate(nodes_generation(nnodes))
          allocate(nodes_area(nnodes))
          allocate(corner_nodes_generation(ncnodes))
          allocate(corner_nodes_area(ncnodes))
          allocate(corner_nodes_faces_gen(ncnodes, 2 * ndims))
          allocate(corner_nodes_faces_area(ncnodes, 2 * ndims))
          allocate(nodes_dirs(nnodes, ndims))
          allocate(corner_nodes_dirs(ncnodes, ndims))
          allocate(face_nodes_dirs(nfnodes, ndims))
          allocate(face_nodes_faces_gen(nfnodes, 2 * ndims))
          allocate(face_nodes_faces_area(nfnodes, 2 * ndims))

          ! Build some uninitialized mesh fields
          call construct_CCSF(mesh_class,  int_CCSF_1)
          call construct_CCSF(mesh_class, real_CCSF_1)

          call construct_CCVF(mesh_class,  int_CCVF_1)
          call construct_CCVF(mesh_class, real_CCVF_1)
          call construct_CCVF(mesh_class,  int_CCVF_2, 2*ndims)
          call construct_CCVF(mesh_class, real_CCVF_2, 2*ndims)

          call construct_FCSF(mesh_class,  int_FCSF_1)
          call construct_FCSF(mesh_class, real_FCSF_1)

          call construct_FCDSF(mesh_class,  int_FCDSF_1)
          call construct_FCDSF(mesh_class, real_FCDSF_1)

          call construct_FCVF(mesh_class,  int_FCVF_1)
          call construct_FCVF(mesh_class, real_FCVF_1)
          call construct_FCVF(mesh_class,  int_FCVF_2, 2*ndims)
          call construct_FCVF(mesh_class, real_FCVF_2, 2*ndims)

          call construct_NCSF(mesh_class,  int_NCSF_1)
          call construct_NCSF(mesh_class, real_NCSF_1)
          call construct_NCSF(mesh_class,  int_NCSF_2, ncnodes)
          call construct_NCSF(mesh_class, real_NCSF_2, ncnodes)

          call construct_NCVF(mesh_class,  int_NCVF_1)
          call construct_NCVF(mesh_class, real_NCVF_1)
          call construct_NCVF(mesh_class,  int_NCVF_2, ncnodes)
          call construct_NCVF(mesh_class, real_NCVF_2, ncnodes)
          call construct_NCVF(mesh_class,  int_NCVF_3, ncnodes, 2*ndims)
          call construct_NCVF(mesh_class, real_NCVF_3, ncnodes, 2*ndims)

          ! Test functions that return large arrays
          nodes_vertices = get_vertices(mesh_class)
          corner_nodes_vertices = get_corner_node_vertices(mesh_class)
          centered_nodes_vertices = get_face_centered_node_vertices(mesh_class)
          cell = 1
          do while (cell .le. ncells)
              ! Test cell-centered values
              volume(cell) = get_cell_volume(mesh_class, cell)
              generation(cell) = get_cell_generation(mesh_class, cell)
              cell_nodes_coords = get_cell_vertices(mesh_class, cell)

              ! Test face-dependent cell values
              face = 1
              do while (face .le. 2 * ndims)
                  num_adj_cell_faces(cell, face) =                      & 
                      get_num_adj(mesh_class, cell, face)
                  adj_cells_faces(cell, face) =                         & 
                      get_next_cell(mesh_class, cell, face)
                  cells_faces_area(cell, face) =                        &
                      get_cell_face_area(mesh_class, cell, face)
                  cell_face_nodes =                                     & 
                      get_cell_face_nodes(mesh_class, cell, face)
                  face_nodes_vertices =                                 &
                      get_cell_face_vertices(mesh_class,cell,face)
                  cell_faces_centered_node(face) =                      &
                      get_cell_face_centered_node(mesh_class, cell, face)
                  cell_faces_area(face) = cells_faces_area(cell, face)
                  cell_faces_specific_node(face) =                      &
                      get_cell_face_centered_node(mesh_class, cell, face)

                  call set_CCVF(int_CCVF_2, cell, face, generation(cell))
                  generation(cell) = get_CCVF(int_CCVF_2, cell, face)
                  call set_CCVF(real_CCVF_2, cell, face, volume(cell))
                  volume(cell) = get_CCVF(real_CCVF_2, cell, face)

                  call set_FCSF(int_FCSF_1, cell, face, generation(cell))
                  generation(cell) = get_FCSF(int_FCSF_1, cell, face)
                  call set_FCSF(real_FCSF_1,cell,face,volume(cell))
                  volume(cell) = get_FCSF(real_FCSF_1, cell, face)

                  call set_FCDSF(int_FCDSF_1, cell, face, generation(cell))
                  generation(cell) = get_FCDSF(int_FCDSF_1, cell, face)
                  call set_FCDSF(real_FCDSF_1, cell, face, volume(cell))
                  volume(cell) = get_FCDSF(real_FCDSF_1, cell, face)

                  call set_FCVF(int_FCVF_1, cell, face, dirs)
                  dirs = get_FCVF(int_FCVF_1, cell, face)
                  call set_FCVF(real_FCVF_1,cell,face, mesh_min_coords)
                  mesh_min_coords = get_FCVF(real_FCVF_1, cell, face)

                  call set_FCVF(int_FCVF_2,cell,face,cell_faces_specific_node)
                  cell_faces_specific_node = get_FCVF(int_FCVF_2, cell, face)
                  call set_FCVF(real_FCVF_2,cell,face, cell_faces_area)
                  cell_faces_area = get_FCVF(real_FCVF_2, cell, face)

                  dir = 1
                  do while (dir .le. ndims)
                      call set_FCVF(int_FCVF_1, cell, face, dir,         &
                          generation(cell))
                      generation(cell) = get_FCVF(int_FCVF_1,cell,face,dir)
                      call set_FCVF(real_FCVF_1,cell,face, dir, volume(cell))
                      volume(cell) = get_FCVF(real_FCVF_1, cell, face, dir)

                      dir = dir + 1
                  end do

                  dir = 1
                  do while (dir .le. 2 * ndims)
                      call set_FCVF(int_FCVF_2, cell, face, dir,         &
                          generation(cell))
                      generation(cell) = get_FCVF(int_FCVF_2,cell,face,dir)
                      call set_FCVF(real_FCVF_2,cell,face, dir, volume(cell))
                      volume(cell) = get_FCVF(real_FCVF_2, cell, face, dir)

                      dir = dir + 1
                  end do

                  face = face + 1

              end do

              ! Test functions that treat all of a single cell's faces
              call set_CCVF(int_CCVF_2, cell, cell_faces_centered_node)
              cell_faces_centered_node = get_CCVF(int_CCVF_2, cell)
              call set_CCVF(real_CCVF_2, cell, cell_faces_area)
              cell_faces_area = get_CCVF(real_CCVF_2, cell)

              call set_FCSF(int_FCSF_1, cell, cell_faces_centered_node)
              cell_faces_centered_node = get_FCSF(int_FCSF_1, cell)
              call set_FCSF(real_FCSF_1, cell, cell_faces_area)
              cell_faces_area = get_FCSF(real_FCSF_1, cell)

              call set_FCDSF(int_FCDSF_1,cell,cell_faces_centered_node)
              cell_faces_centered_node = get_FCDSF(int_FCDSF_1, cell)
              call set_FCDSF(real_FCDSF_1,cell,cell_faces_area)
              cell_faces_area = get_FCDSF(real_FCDSF_1, cell)

              ! Test direction-dependent cell values
              dir = 1
              do while (dir .le. ndims)
                  if (cell .eq. 1) then
                      mesh_min_coords(dir) =                            &
                          get_mesh_min_coord(mesh_class, dir)
                      mesh_max_coords(dir) =                            &
                          get_mesh_max_coord(mesh_class, dir)
                  endif

                  cells_min_coords(cell, dir) =                         &
                      get_cell_min_coord(mesh_class, cell, dir)
                  cells_mid_coords(cell, dir) =                         &
                      get_cell_mid_coord(mesh_class, cell, dir)
                  cells_max_coords(cell, dir) =                         &
                      get_cell_max_coord(mesh_class, cell, dir)
                  cells_dirs_width(cell, dir) =                         &
                      get_cell_width(mesh_class, cell, dir)
                  cells_dirs(cell, dir) = dir

                  call set_CCVF(int_CCVF_1, cell, dir,                  &
                                cells_dirs(cell, dir))
                  cells_dirs(cell, dir) = get_CCVF(int_CCVF_1, cell, dir)
                  call set_CCVF(real_CCVF_1, cell, dir,                 &
                                cells_dirs_width(cell, dir))
                  cells_dirs_width(cell, dir) = get_CCVF(real_CCVF_1,   &
                                                         cell, dir)

                  dirs(dir) = dir

                  dir = dir + 1
              end do

              ! Test node-dependent cell values
              cell_nodes = get_cell_nodes(mesh_class, cell)
              cell_corner_nodes = get_cell_corner_nodes(mesh_class, cell)
              cell_face_centered_nodes =                                &
                  get_cell_face_centered_nodes(mesh_class, cell)
              node = 1
              do while (node .le. (2**ndims + 2*ndims))
                  cells_nodes(cell, node) =                             &
                      get_cell_node(mesh_class, cell, node)
                  node = node + 1
              end do
              node = 1
              do while (node .le. (2**ndims + 2*ndims))
                  cell_nodes_vertices =                                 &
                      get_node_vertices(mesh_class, cell_nodes(node))
                  node = node + 1
              end do

              ! Test mesh field cell-dependent assignment and query operators
              call set_CCSF(int_CCSF_1, cell, generation(cell))
              generation(cell) = get_CCSF(int_CCSF_1, cell)
              call set_CCSF(real_CCSF_1, cell, volume(cell))
              volume(cell) = get_CCSF(real_CCSF_1, cell)

              ! It takes too long to do this for all of the cells just for 
              ! testing
              if (cell .eq. 1) then
                  cell = 0
              endif
              cell = cell + get_num_cells(mesh_class)/4
          end do

          ! Test node-centered functions
          node = 1
          cell = 1
          do while (node .le. nnodes)
              nodes_generation(node) = generation(cell)
              nodes_area(node) = cell_faces_area(cell)

              call set_NCSF(int_NCSF_1, node, generation(cell))
              generation(cell) = get_NCSF(int_NCSF_1, node)
              call set_NCSF(real_NCSF_1, node, nodes_area(node))
              nodes_area(node) = get_NCSF(real_NCSF_1, node)

              call set_NCVF(int_NCVF_1, node, dirs)
              dirs = get_NCVF(int_NCVF_1, node)
              call set_NCVF(real_NCVF_1, node, mesh_min_coords)
              mesh_min_coords = get_NCVF(real_NCVF_1, node)

              dir = 1
              do while (dir .le. ndims)                  
                  call set_NCVF(int_NCVF_1, node, dir, dirs(dir))
                  dirs(dir) = get_NCVF(int_NCVF_1, node, dir)
                  call set_NCVF(real_NCVF_1, node, dir, mesh_min_coords(dir))
                  mesh_min_coords(dir) = get_NCVF(real_NCVF_1, node, dir)

                  nodes_dirs(node, dir) = dir
                  dir = dir + 1
              end do

              if (node .le. ncnodes) then
                  corner_nodes_generation(node) = generation(cell)
                  corner_nodes_area(node) = cell_faces_area(cell)

                  call set_NCSF(int_NCSF_2, node,                       &
                                corner_nodes_generation(node))
                  corner_nodes_generation(node) = get_NCSF(int_NCSF_2, node)
                  call set_NCSF(real_NCSF_2, node, corner_nodes_area(node))
                  corner_nodes_area(node) = get_NCSF(real_NCSF_2, node)

                  call set_NCVF(int_NCVF_2, node, dirs)
                  dirs = get_NCVF(int_NCVF_2, node)
                  call set_NCVF(real_NCVF_2, node, mesh_min_coords)
                  mesh_min_coords = get_NCVF(real_NCVF_2, node)

                  dir = 1
                  do while (dir .le. ndims)
                      call set_NCVF(int_NCVF_2, node, dir, dirs(dir))
                      dirs(dir) = get_NCVF(int_NCVF_2, node, dir)
                      call set_NCVF(real_NCVF_2,node,dir,mesh_min_coords(dir))
                      mesh_min_coords(dir) = get_NCVF(real_NCVF_2, node, dir)
                      corner_nodes_dirs(node, dir) = dir

                      dir = dir + 1
                  end do

                  call set_NCVF(int_NCVF_3, node, cell_face_centered_nodes)
                  cell_face_centered_nodes = get_NCVF(int_NCVF_3, node)
                  call set_NCVF(real_NCVF_3, node, cell_faces_area)
                  cell_faces_area = get_NCVF(real_NCVF_3, node)
                  
                  face = 1
                  do while (face .le. 2*ndims)
                      call set_NCVF(int_NCVF_3, node, face,             &
                                    cell_face_centered_nodes(face))
                      cell_face_centered_nodes(face) =                  &
                          get_NCVF(int_NCVF_3, node, face)
                      call set_NCVF(real_NCVF_3, node, face,            &
                                    cell_faces_area(face))
                      cell_faces_area(face) = get_NCVF(real_NCVF_3,node,face)

                      face = face + 1
                  end do
              endif

              ! It takes too long to do this for all of the nodes just for 
              ! testing
              if (node .eq. 1) then
                  node = 0
              endif

              if (node.lt. ncnodes) then
                  node = node + get_num_corner_nodes(mesh_class)/4
              else
                  node = node + get_num_face_nodes(mesh_class)/4
              end if

          end do

          node = 1
          cell = 1
          face = 1
          do while (cell .le. ncells)
              do while (node .le. 2 * ndims)
                  dis_face_generation(cell, node) = generation(cell)
                  dis_face_nodes_area(cell, node) = cells_faces_area(cell,face)
                  node = node + 1
              end do
              cell = cell + 1
          end do

          ! Build some initialized mesh fields
          call construct_CCSF(mesh_class,   int_CCSF_2, generation)
          call construct_CCSF(mesh_class,  real_CCSF_2, volume)

          call construct_CCVF(mesh_class,   int_CCVF_3, cells_dirs)
          call construct_CCVF(mesh_class,  real_CCVF_3, cells_dirs_width)
          call construct_CCVF(mesh_class,   int_CCVF_4, 2 * ndims,&
                              dis_face_generation)
          call construct_CCVF(mesh_class,  real_CCVF_4, 2 * ndims,&
                              dis_face_nodes_area)

          call construct_FCSF(mesh_class,  int_FCSF_2, face_generation)
          call construct_FCSF(mesh_class, real_FCSF_2, face_nodes_area)

          call construct_FCDSF(mesh_class,  int_FCDSF_2,          &
                               dis_face_generation)
          call construct_FCDSF(mesh_class, real_FCDSF_2,          &
                               dis_face_nodes_area)

          call construct_FCVF(mesh_class,  int_FCVF_3, face_nodes_dirs)
          call construct_FCVF(mesh_class, real_FCVF_3, centered_nodes_vertices)
          call construct_FCVF(mesh_class,  int_FCVF_4, face_nodes_faces_gen)
          call construct_FCVF(mesh_class, real_FCVF_4, face_nodes_faces_area)

          call construct_NCSF(mesh_class,  int_NCSF_3, nodes_generation)
          call construct_NCSF(mesh_class, real_NCSF_3, nodes_area)
          call construct_NCSF(mesh_class,  int_NCSF_4,            &
                              corner_nodes_generation)
          call construct_NCSF(mesh_class,  real_NCSF_4,           &
                              corner_nodes_area)

          call construct_NCVF(mesh_class,  int_NCVF_4, nodes_dirs)
          call construct_NCVF(mesh_class, real_NCVF_4, nodes_vertices)
          call construct_NCVF(mesh_class,  int_NCVF_5,            &
                              corner_nodes_dirs)
          call construct_NCVF(mesh_class,  real_NCVF_5,           &
                              corner_nodes_vertices)
          call construct_NCVF(mesh_class,  int_NCVF_6,corner_nodes_faces_gen)
          call construct_NCVF(mesh_class, real_NCVF_6,corner_nodes_faces_area)

          ! Test mesh field assignment and query operators
          call set_CCSF(int_CCSF_1, generation)
          generation = get_CCSF(int_CCSF_1)
          call set_CCSF(real_CCSF_1, volume)
          volume = get_CCSF(real_CCSF_1)

          call set_CCVF(int_CCVF_4, dis_face_generation)
          dis_face_generation = get_CCVF(int_CCVF_4)
          call set_CCVF(real_CCVF_4, dis_face_nodes_area)
          dis_face_nodes_area = get_CCVF(real_CCVF_4)

          call set_FCSF(int_FCSF_1, face_generation)
          face_generation = get_FCSF(int_FCSF_1)
          call set_FCSF(real_FCSF_1, face_nodes_area)
          face_nodes_area  = get_FCSF(real_FCSF_1)

          call set_FCDSF(int_FCDSF_1, dis_face_generation)
          dis_face_generation = get_FCDSF(int_FCDSF_1)
          call set_FCDSF(real_FCDSF_1, dis_face_nodes_area)
          dis_face_nodes_area = get_FCDSF(real_FCDSF_1)

          call set_FCVF( int_FCVF_3, face_nodes_dirs)
          face_nodes_dirs = get_FCVF( int_FCVF_3)
          call set_FCVF(real_FCVF_3, centered_nodes_vertices)
          centered_nodes_vertices = get_FCVF(real_FCVF_3)
          call set_FCVF( int_FCVF_4, face_nodes_faces_gen)
          face_nodes_faces_gen = get_FCVF( int_FCVF_4)
          call set_FCVF(real_FCVF_4, face_nodes_faces_area)
          face_nodes_faces_area = get_FCVF(real_FCVF_4)

          call set_NCSF(int_NCSF_3, nodes_generation)
          nodes_generation = get_NCSF(int_NCSF_3)
          call set_NCSF(real_NCSF_3, nodes_area)
          nodes_area = get_NCSF(real_NCSF_3)
          call set_NCSF(int_NCSF_4, corner_nodes_generation)
          corner_nodes_generation = get_NCSF(int_NCSF_4)
          call set_NCSF(real_NCSF_4, corner_nodes_area)
          corner_nodes_area = get_NCSF(real_NCSF_4)

          call set_NCVF( int_NCVF_4, nodes_dirs)
          nodes_dirs = get_NCVF( int_NCVF_4)
          call set_NCVF(real_NCVF_4, nodes_vertices)
          nodes_vertices = get_NCVF(real_NCVF_4)
          call set_NCVF(int_NCVF_5, corner_nodes_dirs)
          corner_nodes_dirs = get_NCVF(int_NCVF_5)
          call set_NCVF(real_NCVF_5,corner_nodes_vertices)
          corner_nodes_vertices = get_NCVF(real_NCVF_5)
          call set_NCVF(int_NCVF_6,corner_nodes_faces_gen)
          corner_nodes_faces_gen = get_NCVF(int_NCVF_6)
          call set_NCVF(real_NCVF_6,corner_nodes_faces_area)
          corner_nodes_faces_area = get_NCVF(real_NCVF_6)

          ! Deallocate memory for test variable arrays
          ! Cell-centered values
          deallocate(volume)
          deallocate(generation)
          ! Face-dependent cell values
          deallocate(num_adj_cell_faces)
          deallocate(adj_cells_faces)
          deallocate(cells_faces_area)
          deallocate(cell_faces_centered_node)
          ! Face-dependent values
          deallocate(cell_faces_area)
          ! Direction-dependent values
          deallocate(mesh_min_coords)
          deallocate(mesh_max_coords)
          deallocate(dirs)
          ! Direction-dependent cell values
          deallocate(cells_min_coords)
          deallocate(cells_mid_coords)
          deallocate(cells_max_coords)
          deallocate(cells_dirs_width)
          deallocate(cells_dirs)
          ! Node-dependent cell values
          deallocate(cell_face_nodes)
          deallocate(cell_faces_specific_node)
          deallocate(cell_nodes)
          deallocate(cell_corner_nodes)
          deallocate(cell_face_centered_nodes)
          deallocate(cells_nodes)
          ! Node-centered values
          deallocate(nodes_vertices)
          deallocate(corner_nodes_vertices)
          deallocate(centered_nodes_vertices)
          deallocate(face_nodes_vertices)
          deallocate(cell_nodes_coords)
          deallocate(cell_nodes_vertices)
          deallocate(face_generation)
          deallocate(face_nodes_area)
          deallocate(dis_face_generation)
          deallocate(dis_face_nodes_area)
          deallocate(nodes_generation)
          deallocate(nodes_area)
          deallocate(corner_nodes_generation)
          deallocate(corner_nodes_area)
          deallocate(corner_nodes_faces_gen)
          deallocate(corner_nodes_faces_area)
          deallocate(nodes_dirs)
          deallocate(corner_nodes_dirs)
          deallocate(face_nodes_dirs)
          deallocate(face_nodes_faces_gen)
          deallocate(face_nodes_faces_area)

!===========================================================================
! Get rid of the CAR_CU_Mesh class object and the associated test fields. 
! We are done.
!===========================================================================

          call destruct_CCSF(int_CCSF_1)
          call destruct_CCSF(int_CCSF_2)
          call destruct_CCSF(real_CCSF_1)
          call destruct_CCSF(real_CCSF_2)

          call destruct_CCVF(int_CCVF_1)
          call destruct_CCVF(int_CCVF_2)
          call destruct_CCVF(int_CCVF_3)
          call destruct_CCVF(int_CCVF_4)
          call destruct_CCVF(real_CCVF_1)
          call destruct_CCVF(real_CCVF_2)
          call destruct_CCVF(real_CCVF_3)
          call destruct_CCVF(real_CCVF_4)

          call destruct_FCSF(int_FCSF_1)
          call destruct_FCSF(int_FCSF_2)
          call destruct_FCSF(real_FCSF_1)
          call destruct_FCSF(real_FCSF_2)

          call destruct_FCDSF(int_FCDSF_1)
          call destruct_FCDSF(int_FCDSF_2)
          call destruct_FCDSF(real_FCDSF_1)
          call destruct_FCDSF(real_FCDSF_2)

          call destruct_FCVF(int_FCVF_1)
          call destruct_FCVF(int_FCVF_2)
          call destruct_FCVF(int_FCVF_3)
          call destruct_FCVF(int_FCVF_4)
          call destruct_FCVF(real_FCVF_1)
          call destruct_FCVF(real_FCVF_2)
          call destruct_FCVF(real_FCVF_3)
          call destruct_FCVF(real_FCVF_4)

          call destruct_NCSF(int_NCSF_1)
          call destruct_NCSF(int_NCSF_2)
          call destruct_NCSF(int_NCSF_3)
          call destruct_NCSF(int_NCSF_4)
          call destruct_NCSF(real_NCSF_1)
          call destruct_NCSF(real_NCSF_2)
          call destruct_NCSF(real_NCSF_3)
          call destruct_NCSF(real_NCSF_4)

          call destruct_NCVF(int_NCVF_1)
          call destruct_NCVF(int_NCVF_2)
          call destruct_NCVF(int_NCVF_3)
          call destruct_NCVF(int_NCVF_4)
          call destruct_NCVF(int_NCVF_5)
          call destruct_NCVF(int_NCVF_6)
          call destruct_NCVF(real_NCVF_1)
          call destruct_NCVF(real_NCVF_2)
          call destruct_NCVF(real_NCVF_3)
          call destruct_NCVF(real_NCVF_4)
          call destruct_NCVF(real_NCVF_5)
          call destruct_NCVF(real_NCVF_6)

          call destruct_Mesh(mesh_class)

      end program alpha

