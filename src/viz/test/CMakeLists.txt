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
file( GLOB python_scripts python/*.py )
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

add_scalar_tests(
   SOURCES "${test_sources}"
   DEPS    "Lib_traits;Lib_c4;Lib_dsxx;${MPI_LIBRARIES};${INTEL_LIBM}"
)

# ---------------------------------------------------------------------------- #
# Copy input files to test directory (binary_dir)
# ---------------------------------------------------------------------------- #

provide_aux_files(
   FILES "${mesh_input_files}"
   SRC_EXT ".in"
)
provide_aux_files( FILES "${python_scripts}" )
provide_aux_files( FILES "${bench_files}" )


