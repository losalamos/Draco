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
set( test_deps
   Lib_dsxx
   ${MPI_LIBRARIES} 
   ${PAPI_LIBRARY}
   ${omplib}
)


# This test has special arguments and special pass/fail criteria.
set( special_test ${PROJECT_SOURCE_DIR}/phw.cc )
add_parallel_tests( 
   SOURCES   "${special_test}" 
   PE_LIST   "1;2;4" 
   DEPS      "${test_deps}"
   TEST_ARGS "hello"
   PASS_REGEX "Hello"   
   FAIL_REGEX "cruel"
   )
list( REMOVE_ITEM test_sources ${special_test} )

# ApplicationUnitTest uses a paradigm that is not compatbilty with
# Catamount systems.  We cannot run aprun under a system call.
set( special_test ${PROJECT_SOURCE_DIR}/tstApplicationUnitTest.cc )
if( NOT "${CMAKE_SYSTEM_NAME}" MATCHES "Catamount" )
   # This test also has special arguments and spawns mpirun indirectly so
   # we register it as a scalar test.
   if( "${DRACO_C4}" MATCHES "MPI" )
      set( application_unit_test_args "--np 1;--np 2;--np 4" )
   else()
      set( application_unit_test_args "--np scalar" )
   endif()
   add_scalar_tests(
      SOURCES   "${special_test}" 
      DEPS      "${test_deps}"
      TEST_ARGS "${application_unit_test_args}"
      )
endif()
list( REMOVE_ITEM test_sources ${special_test} )

# This tests runs with extra arguments 
set( special_test ${PROJECT_SOURCE_DIR}/tstParallelUnitTest.cc )
add_parallel_tests(
   SOURCES   "${special_test}" 
   DEPS      "${test_deps}"
   PE_LIST   "1;2;4"
   TEST_ARGS "--version;--dummy; "
   PASS_REGEX "tstParallelUnitTest: version"
   )
list( REMOVE_ITEM test_sources ${special_test} )

# The rest of the tests

add_parallel_tests( 
   SOURCES   "${test_sources}" 
   PE_LIST   "1;2;4" 
   DEPS      "${test_deps}"
   )
   
# tstOMP should be run alone (not in paralle with other tests) because
# it tries to use all available cores as OMP threads.
set_tests_properties( c4_tstOMP_1 c4_tstOMP_2 c4_tstOMP_4 
   PROPERTIES APPEND
      RUN_SERIAL
   )

 # RUN_SERIAL
 #       Do not run this test in parallel with any other test.
 #       Use this option in conjunction with the ctest_test PARALLEL_LEVEL
 #       option to specify that this test should not be run in parallel with
 #       any other tests.
