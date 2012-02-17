# ---------------------------------------------------------------------------- #
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
# ---------------------------------------------------------------------------- #
# $Id$
# ---------------------------------------------------------------------------- #

project( viz_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_lib_sources 
   ${PROJECT_SOURCE_DIR}/viz_test.cc 
   ${PROJECT_SOURCE_DIR}/viz_test.hh
)
file( GLOB test_sources *.cc )
list( REMOVE_ITEM test_sources ${test_lib_sources} )
file( GLOB mesh_input_files *.in )
file( GLOB bench_files bench/*.0001 )

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories( 
   ${PROJECT_SOURCE_DIR}      # headers for tests
   ${PROJECT_SOURCE_DIR}/..   # headers for package
)

# ---------------------------------------------------------------------------- #
# Build library for test directory
# ---------------------------------------------------------------------------- #

add_component_library( Lib_viz_test ${PROJECT_NAME} "${test_lib_sources}" )
add_dependencies( Lib_viz_test Lib_viz )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

# This test requires numdiff, so don't run it if it's not available.
find_program( NUMDIFF
              numdiff
              HINTS ${VENDOR_DIR}
            )

#if( NOT "${CMAKE_CXX_COMPILER}" MATCHES "[sp]pu-g[+][+]" )
if ( NOT ${NUMDIFF} MATCHES NOTFOUND )
   add_scalar_tests(
      SOURCES "${test_sources}"
      DEPS    "Lib_dsxx" 
    # "Lib_c4;${MPI_LIBRARIES}"
      )
endif()

# ---------------------------------------------------------------------------- #
# Copy input files to test directory (binary_dir)
# ---------------------------------------------------------------------------- #

provide_aux_files(
   FILES "${mesh_input_files}"
   SRC_EXT ".in"
)
provide_aux_files( FILES "${python_scripts}" )
provide_aux_files( FILES "${bench_files}" )


##---------------------------------------------------------------------------##
## Extra 'make clean' files
##---------------------------------------------------------------------------##
set( extra_clean_files
   Densities.diff
   Pressure.diff
   Temperatures.diff
   Velocity.diff
   ensight_stream.out
   geo.diff
   part_testproblem_binary_ensight
   part_testproblem_ensight
   testproblem_binary_ensight
   testproblem_ensight
)
set_directory_properties(
   PROPERTIES
   ADDITIONAL_MAKE_CLEAN_FILES "${extra_clean_files}" )


