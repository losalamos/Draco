# ---------------------------------------------------------------------------- #
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
# ---------------------------------------------------------------------------- #
# $Id$
# ---------------------------------------------------------------------------- #

project( meshReaders_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_lib_sources 
   ${PROJECT_SOURCE_DIR}/meshReaders_test.cc 
   ${PROJECT_SOURCE_DIR}/meshReaders_test.hh
)
file( GLOB test_sources *.cc )
file( GLOB mesh_input_files *.mesh.in )
list( REMOVE_ITEM test_sources ${test_lib_sources} )

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

add_component_library( Lib_meshReaders_test ${PROJECT_NAME} "${test_lib_sources}" )
add_dependencies( Lib_meshReaders_test Lib_meshReaders )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_scalar_tests(
   SOURCES "${test_sources}"
   DEPS    "Lib_dsxx;${INTEL_LIBM}" )

# ---------------------------------------------------------------------------- #
# Copy mesh files to test directory (binary_dir)
# ---------------------------------------------------------------------------- #

provide_aux_files(
   FILES    "${mesh_input_files}"
   SRC_EXT  ".mesh.in"
   DEST_EXT ".mesh" )


