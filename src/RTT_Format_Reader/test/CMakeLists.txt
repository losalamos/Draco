#=============================================================================
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
#=============================================================================
# $Id$
#=============================================================================

project( RTT_Format_Reader_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_lib_sources 
	${PROJECT_SOURCE_DIR}/RTT_Format_Reader_test.cc 
	${PROJECT_SOURCE_DIR}/RTT_Format_Reader_test.hh)
file( GLOB test_sources *.cc )
list( REMOVE_ITEM test_sources ${test_lib_sources} )
file( GLOB mesh_input_files *.mesh.in )

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories( 
   ${PROJECT_SOURCE_DIR}      # headers for tests
   ${PROJECT_SOURCE_DIR}/..   # headers for package
   ${PROJECT_BINARY_DIR}/..   # config.h
   ${MPI_INCLUDE_PATH} )

# ---------------------------------------------------------------------------- #
# Build library for test directory
# ---------------------------------------------------------------------------- #

add_component_library( Lib_RTT_Format_Reader_test 
   RTT_Format_Reader_test 
   ${test_lib_sources} )
add_dependencies( Lib_RTT_Format_Reader_test
   Lib_RTT_Format_Reader 
   Lib_RTT_Format_Reader
   Lib_mesh_element
   Lib_c4
   Lib_dsxx)

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

set(test_deps
   Lib_RTT_Format_Reader_test 
   Lib_RTT_Format_Reader
   Lib_mesh_element
   Lib_c4 
   Lib_dsxx 
   ${MPI_LIBRARIES} )

add_scalar_tests(
   SOURCES "${test_sources}"
   DEPS    "${test_deps}" )

# ---------------------------------------------------------------------------- #
# Copy input files to test directory (binary_dir)
# ---------------------------------------------------------------------------- #

provide_aux_files(
   FILES    "${mesh_input_files}"
   SRC_EXT  ".mesh.in"
   DEST_EXT ".mesh" )

