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
!
! Shadow Interface modules
!
!===========================================================================

          USE CAR_CU_Interface_Class
          USE CAR_CU_RTT_Format_Class
          USE CAR_CU_Builder_Class
          USE CAR_CU_Mesh_Class

          implicit none
          integer narg, iargc, fnlgth
          type(CAR_CU_Interface) :: interface_class
          type(CAR_CU_RTT_Format) :: rtt_format_class
          type(CAR_CU_Builder) :: builder_class
          type(CAR_CU_Mesh) :: mesh_class

!===========================================================================
! Define variables just needed for testing the shadow interface functions
!===========================================================================
          integer ndims, dir, ncells, nnodes, ncnodes, nfnodes, cell,   &
              face, node_index
          integer, dimension (:,:), allocatable :: num_adj, adj_cell,   &
              face_area, cell_specific_nodes
          real, dimension (:,:), allocatable    :: cell_min_val,        &
              cell_mid_val,cell_max_val          
          integer, dimension (:), allocatable   :: generation,          &
              cell_nodes, cell_corner_nodes, cell_face_cen_nodes,       &
              cell_face_nodes, cell_face_specific_nodes

          real, dimension (:), allocatable     :: volume, mesh_min_val, &
              mesh_max_val

!===========================================================================
! Input the command line arguments - input file name followed by anything to
! activate the verbose switch (typically v or the word verbose)
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
! that are no longer needed. 
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
          ! Direction-dependent values
          allocate(mesh_min_val(ndims))
          allocate(mesh_max_val(ndims))
          ! Direction-dependent cell values
          allocate(cell_min_val(ncells, ndims))
          allocate(cell_mid_val(ncells, ndims))
          allocate(cell_max_val(ncells, ndims))
          ! Node-dependent cell values
          allocate(cell_specific_nodes(ncells, (2**ndims + 2*ndims)))
          allocate(cell_nodes(2**ndims + 2*ndims))
          allocate(cell_corner_nodes(2**ndims))
          allocate(cell_face_cen_nodes(2 * ndims))
          allocate(cell_face_nodes(2 * (ndims -1)))
          allocate(cell_face_specific_nodes(2 * (ndims -1)))

          cell = 1
          do while (cell .le. ncells)
              ! Test cell-centered values
              volume(cell) = get_cell_volume(mesh_class, cell)
              generation(cell) = get_cell_generation(mesh_class, cell)

              ! Test node-dependent cell values
              cell_nodes = get_cell_nodes(mesh_class, cell)
              cell_corner_nodes = get_cell_corner_nodes(mesh_class, cell)
              cell_face_cen_nodes =                                     &
                  get_cell_face_centered_nodes(mesh_class, cell)
              node_index = 1
              do while (node_index .le. (2**ndims + 2*ndims))
                 cell_specific_nodes(cell, node_index) =                &
                     get_cell_node(mesh_class, cell, node_index)

                 node_index = node_index + 1
              end do

              ! Test face-dependent cell values
              face = 1
              do while (face .le. 2 * ndims)
                  num_adj(cell, face) = get_num_adj(mesh_class, cell, face)
                  adj_cell(cell, face) = get_next_cell(mesh_class, cell, face)
                  face_area(cell, face) =                               &
                      get_cell_face_area(mesh_class, cell, face)
                  cell_face_nodes = get_cell_face_nodes(mesh_class, cell, face)
                  cell_face_specific_nodes(face) =                      &
                      get_cell_face_centered_node(mesh_class, cell, face)

                  face = face + 1
              end do

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

                  dir = dir + 1
              end do
              cell = cell + 1
          end do

          ! Deallocate memory for test variable arrays
          ! Cell-centered values
          deallocate(volume)
          deallocate(generation)
          ! Face-dependent cell values
          deallocate(num_adj)
          deallocate(adj_cell)
          deallocate(face_area)
          ! Direction-dependent values
          deallocate(mesh_min_val)
          deallocate(mesh_max_val)
          ! Direction-dependent cell values
          deallocate(cell_min_val)
          deallocate(cell_mid_val)
          deallocate(cell_max_val)
          ! Node-dependent cell values
          deallocate(cell_specific_nodes)
          deallocate(cell_nodes)
          deallocate(cell_corner_nodes)
          deallocate(cell_face_cen_nodes)
          deallocate(cell_face_nodes)
          deallocate(cell_face_specific_nodes)


!===========================================================================
! Get rid of the CAR_CU_Mesh class object - we are done.
!===========================================================================

          call destruct_Mesh_Class(mesh_class)

      end program alpha

