#=============================================================================
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
#=============================================================================
# $Id$
#=============================================================================

project( traits_test CXX )

file( GLOB test_sources *.cc )

# Find local headers
include_directories( 
   ${PROJECT_SOURCE_DIR}      # headers for tests
   ${PROJECT_SOURCE_DIR}/..   # headers for package
   ${dsxx_BINARY_DIR}         # headers for ds++
   ${draco_src_dir_SOURCE_DIR}
   ${MPI_INCLUDE_PATH} )

add_scalar_tests(
   SOURCES "${test_sources}"
   DEPS    "Lib_dsxx"
   )

