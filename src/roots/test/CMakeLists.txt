# ---------------------------------------------------------------------------- #
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
# ---------------------------------------------------------------------------- #
# $Id$
# ---------------------------------------------------------------------------- #

project( roots_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_lib_sources 
   ${PROJECT_SOURCE_DIR}/roots_test.cc 
   ${PROJECT_SOURCE_DIR}/roots_test.hh
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
   ${c4_BINARY_DIR}           # c4/config.h
   ${MPI_INCLUDE_PATH}
)

# ---------------------------------------------------------------------------- #
# Build library for test directory
# ---------------------------------------------------------------------------- #

add_component_library( Lib_roots_test ${PROJECT_NAME} "${test_lib_sources}" )
add_dependencies( Lib_roots_test Lib_roots )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_scalar_tests(
   SOURCES "${test_sources}"
   DEPS "Lib_linear;Lib_c4;Lib_dsxx;${MPI_LIBRARIES};${INTEL_LIBM}" )
