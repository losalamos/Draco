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

          type(integer_CCSF)      :: iCCSF_class, idCCSF_class
          type(real_CCSF)         :: rCCSF_class, rdCCSF_class
          type(integer_FCSF)      :: iFCSF_class, idFCSF_class
          type(real_FCSF)         :: rFCSF_class, rdFCSF_class
          integer ndims, dir, ncells, nnodes, ncnodes, nfnodes, cell,   &
              face, node
          integer, dimension (:,:), allocatable :: num_adj, adj_cell,   &
              face_area, cell_specific_nodes
          real*8, dimension (:,:), allocatable  :: cell_min_val,        &
              cell_mid_val, cell_max_val, cell_width, vertices,         &
              corner_node_vertices, face_centered_node_vertices,        &
              cell_vertices, face_vertices
          integer, dimension (:), allocatable   :: generation,          &
              cell_nodes, cell_corner_nodes, cell_face_cen_nodes,       &
              cell_face_nodes, cell_face_specific_nodes, int_face_nodes,&
              cell_face_node_nums
          real*8, dimension (:), allocatable    :: volume,mesh_min_val, &
              mesh_max_val, node_vertices, real_face_nodes,             &
              cell_face_areas

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

          call construct_Interface_Class(interface_class)
          rtt_format_class%this = interface_class%rtt_format

!===========================================================================
! Create a C++ CAR_CU_Builder class object. This also constructs the C++ 
! Coord_sys, Layout, and CAR_CU_Mesh class objects. The addresses of both 
! the new CAR_CU_Builder and CAR_CU_Mesh class objects are set.
!===========================================================================

          call construct_Builder_Class(builder_class, interface_class)
          mesh_class%this = builder_class%mesh

!===========================================================================
! Get rid of RTT_Format, CAR_CU_Interface, and CAR_CU_Builder class objects 
! since they are no longer needed. 
!===========================================================================

          call destruct_RTT_Format_Class(rtt_format_class)
          call destruct_Interface_Class(interface_class)
          call destruct_Builder_Class(builder_class)

!===========================================================================
! Test shadowed mesh accessor functions
!===========================================================================

          ! Test scalar values
          ndims = get_dimension(mesh_class)
          ncells = get_num_cells(mesh_class)
          nnodes = get_num_nodes(mesh_class)
          ncnodes = get_num_corner_nodes(mesh_class)
          nfnodes = get_num_face_nodes(mesh_class)

          ! Allocate memory for test variable arrays
          ! Cell-centered values
          allocate(volume(ncells))
          allocate(generation(ncells))
          ! Face-dependent cell values
          allocate(num_adj(ncells, 2 * ndims))
          allocate(adj_cell(ncells, 2 * ndims))
          allocate(face_area(ncells, 2 * ndims))
          allocate(cell_face_node_nums( 2 * ndims))
          allocate(cell_face_areas( 2 * ndims))
          ! Direction-dependent values
          allocate(mesh_min_val(ndims))
          allocate(mesh_max_val(ndims))
          ! Direction-dependent cell values
          allocate(cell_min_val(ncells, ndims))
          allocate(cell_mid_val(ncells, ndims))
          allocate(cell_max_val(ncells, ndims))
          allocate(cell_width(ncells, ndims))
          ! Node-dependent cell values
          allocate(cell_specific_nodes(ncells, (2**ndims + 2*ndims)))
          allocate(cell_nodes(2**ndims + 2*ndims))
          allocate(cell_corner_nodes(2**ndims))
          allocate(cell_face_cen_nodes(2 * ndims))
          allocate(cell_face_nodes(2 * (ndims -1)))
          allocate(cell_face_specific_nodes(2 * ndims))
          ! Node-centered values
          allocate(vertices(nnodes, ndims))
          allocate(corner_node_vertices(ncnodes, ndims))
          allocate(face_centered_node_vertices(nfnodes, ndims))
          allocate(cell_vertices((2**ndims), ndims))
          allocate(face_vertices((2*(ndims - 1)), ndims))
          allocate(node_vertices(ndims))
          allocate(int_face_nodes(nfnodes))
          allocate(real_face_nodes(nfnodes))

          ! Build some uninitialized mesh fields
          call construct_integer_CCSF_Class(mesh_class, iCCSF_Class)
          call construct_real_CCSF_Class(mesh_class, rCCSF_Class)
          call construct_integer_FCSF_Class(mesh_class, iFCSF_Class)
          call construct_real_FCSF_Class(mesh_class, rFCSF_Class)

          ! Test functions that return large arrays
          vertices = get_vertices(mesh_class)
          corner_node_vertices = get_corner_node_vertices(mesh_class)
          face_centered_node_vertices =                                 &
                                   get_face_centered_node_vertices(mesh_class)
          cell = 1
          do while (cell .le. ncells)
              ! Test cell-centered values
              volume(cell) = get_cell_volume(mesh_class, cell)
              generation(cell) = get_cell_generation(mesh_class, cell)
              cell_vertices = get_cell_vertices(mesh_class, cell)

              ! Test face-dependent cell values
              face = 1
              do while (face .le. 2 * ndims)
                  num_adj(cell, face) = get_num_adj(mesh_class, cell, face)
                  adj_cell(cell, face) = get_next_cell(mesh_class, cell, face)
                  face_area(cell, face) =                               &
                      get_cell_face_area(mesh_class, cell, face)
                  cell_face_nodes = get_cell_face_nodes(mesh_class, cell, face)
                  face_vertices = get_cell_face_vertices(mesh_class,cell,face)

                  cell_face_node_nums(face) =                           &
                      get_cell_face_centered_node(mesh_class, cell, face)
                  cell_face_areas(face) = face_area(cell, face)
                  cell_face_specific_nodes(face) =                      &
                      get_cell_face_centered_node(mesh_class, cell, face)
                  call set_integer_FCSF_cell_face(iFCSF_Class, cell,    &
                                                  face, generation(cell))
                  generation(cell) =                                    &
                      get_integer_FCSF_cell_face(iFCSF_Class, cell, face)
                  call set_real_FCSF_cell_face(rFCSF_Class, cell, face, &
                      volume(cell))
                  volume(cell) =                                        &
                      get_real_FCSF_cell_face(rFCSF_Class, cell, face)

                  face = face + 1
              end do

              ! Test functions that treat all of a cells faces
              call set_integer_FCSF_cell(iFCSF_Class,cell,cell_face_node_nums)
              cell_face_node_nums = get_integer_FCSF_cell(iFCSF_Class, cell)
              call set_real_FCSF_cell(rFCSF_Class, cell, cell_face_areas)
              cell_face_areas = get_real_FCSF_cell(rFCSF_Class, cell)

              ! Test direction-dependent cell values
              dir = 1
              do while (dir .le. ndims)
                  if (cell .eq. 1) then
                      mesh_min_val(dir) = get_mesh_min_coord(mesh_class, dir)
                      mesh_max_val(dir) = get_mesh_max_coord(mesh_class, dir)
                  endif

                  cell_min_val(cell,dir) =                              &
                      get_cell_min_coord(mesh_class, cell, dir)
                  cell_mid_val(cell,dir) =                              &
                      get_cell_mid_coord(mesh_class, cell, dir)
                  cell_max_val(cell,dir) =                              &
                      get_cell_max_coord(mesh_class, cell, dir)
                  cell_width(cell,dir) = get_cell_width(mesh_class, cell, dir)

                  dir = dir + 1
              end do

              ! Test node-dependent cell values
              cell_nodes = get_cell_nodes(mesh_class, cell)
              cell_corner_nodes = get_cell_corner_nodes(mesh_class, cell)
              cell_face_cen_nodes =                                     &
                  get_cell_face_centered_nodes(mesh_class, cell)
              node = 1
              do while (node .le. (2**ndims + 2*ndims))
                 cell_specific_nodes(cell, node) =                      &
                     get_cell_node(mesh_class, cell, node)

                 node = node + 1
              end do
              node = 1
              do while (node .le. (2**ndims + 2*ndims))
                 node_vertices =                                        &
                     get_node_vertices(mesh_class, cell_nodes(node))

                 node = node + 1
              end do

              ! Test mesh field cell-dependent assignment and query operators
              call set_integer_CCSF_cell(iCCSF_Class, cell, generation(cell))
              generation(cell) = get_integer_CCSF_cell(iCCSF_Class, cell)
              call set_real_CCSF_cell(rCCSF_Class, cell, volume(cell))
              volume(cell) = get_real_CCSF_cell(rCCSF_Class, cell)

              cell = cell + 1
          end do

          ! Build some initialized mesh fields
          call construct_integer_CCSF_Class(mesh_class, idCCSF_Class,   &
                                            generation)
          call construct_real_CCSF_Class(mesh_class, rdCCSF_Class,      &
                                         volume)
          node = 1
          cell = 1
          do while (node .le. nfnodes)
              int_face_nodes(node) = generation(cell)
              real_face_nodes(node) = volume(cell)

              node = node + 1
          end do
          call construct_integer_FCSF_Class(mesh_class, idFCSF_Class,   &
                                            int_face_nodes)
          call construct_real_FCSF_Class(mesh_class, rdFCSF_Class,      &
                                         real_face_nodes)


          ! Test mesh field assignment and query operators
          call set_integer_CCSF(iCCSF_Class, generation)
          generation = get_integer_CCSF(iCCSF_Class)
          call set_real_CCSF(rCCSF_Class, volume)
          volume = get_real_CCSF(rCCSF_Class)

          call set_integer_FCSF(iFCSF_Class, int_face_nodes)
          int_face_nodes = get_integer_FCSF(iFCSF_Class)

          call set_real_FCSF(rFCSF_Class, real_face_nodes)
          real_face_nodes = get_real_FCSF(rFCSF_Class)

          ! Deallocate memory for test variable arrays
          ! Cell-centered values
          deallocate(volume)
          deallocate(generation)
          ! Face-dependent cell values
          deallocate(num_adj)
          deallocate(adj_cell)
          deallocate(face_area)
          deallocate(cell_face_node_nums)
          deallocate(cell_face_areas)
          ! Direction-dependent values
          deallocate(mesh_min_val)
          deallocate(mesh_max_val)
          ! Direction-dependent cell values
          deallocate(cell_min_val)
          deallocate(cell_mid_val)
          deallocate(cell_max_val)
          deallocate(cell_width)
          ! Node-dependent cell values
          deallocate(cell_specific_nodes)
          deallocate(cell_nodes)
          deallocate(cell_corner_nodes)
          deallocate(cell_face_cen_nodes)
          deallocate(cell_face_nodes)
          deallocate(cell_face_specific_nodes)
          ! Node-centered values
          deallocate(vertices)
          deallocate(corner_node_vertices)
          deallocate(face_centered_node_vertices)
          deallocate(cell_vertices)
          deallocate(face_vertices)
          deallocate(node_vertices)
          deallocate(int_face_nodes)
          deallocate(real_face_nodes)


!===========================================================================
! Get rid of the CAR_CU_Mesh class object and the associated fields. 
! We are done.
!===========================================================================

          call destruct_integer_CCSF_Class(iCCSF_Class)
          call destruct_integer_CCSF_Class(idCCSF_Class)
          call destruct_real_CCSF_Class(rCCSF_Class)
          call destruct_real_CCSF_Class(rdCCSF_Class)
          call destruct_Mesh_Class(mesh_class)

      end program alpha

