# ---------------------------------------------------------------------------- #
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
# ---------------------------------------------------------------------------- #
# $Id$
# ---------------------------------------------------------------------------- #

project( c4_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_lib_sources 
   ${PROJECT_SOURCE_DIR}/c4_test.cc 
   ${PROJECT_SOURCE_DIR}/c4_test.hh
)
file( GLOB test_sources *.cc )
list( REMOVE_ITEM test_sources ${test_lib_sources} )

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

add_component_library( Lib_c4_test c4_test "${test_lib_sources}" )
add_dependencies( Lib_c4_test Lib_c4 )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

# This test has special arguments and special pass/fail criteria.
set( special_test ${PROJECT_SOURCE_DIR}/phw.cc )
add_parallel_tests( 
   SOURCES   "${special_test}" 
   PE_LIST   "1;2;4" 
   DEPS      "Lib_dsxx;${MPI_LIBRARIES}"
   TEST_ARGS "hello"
   PASS_REGEX "Hello"   
   FAIL_REGEX "cruel"
   )
list( REMOVE_ITEM test_sources ${special_test} )

# This test also has special arguments and spawns mpirun indirectly so
# we register it as a scalar test.
set( special_test ${PROJECT_SOURCE_DIR}/tstApplicationUnitTest.cc )
add_scalar_tests(
   SOURCES   "${special_test}" 
   DEPS      "Lib_dsxx;${MPI_LIBRARIES}"
   TEST_ARGS "--np 1;--np 2;--np 4"
   )
list( REMOVE_ITEM test_sources ${special_test} )

# This tests runs with extra arguments 
set( special_test ${PROJECT_SOURCE_DIR}/tstParallelUnitTest.cc )
add_parallel_tests(
   SOURCES   "${special_test}" 
   DEPS      "Lib_dsxx;${MPI_LIBRARIES}"
   PE_LIST   "1;2;4"
   TEST_ARGS "--version;--dummy; "
   PASS_REGEX "tstParallelUnitTest: version"
   )
list( REMOVE_ITEM test_sources ${special_test} )

# The rest of the tests
add_parallel_tests( 
   SOURCES   "${test_sources}" 
   PE_LIST   "1;2;4" 
   DEPS      "Lib_dsxx;${MPI_LIBRARIES}"
   )
   

