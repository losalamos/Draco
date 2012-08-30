#-----------------------------*-cmake-*----------------------------------------#
# file   config/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for c4/test.
# note   Copyright (C) 2010-2012, Los Alamos National Security
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#
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

# Bullseye doesn't like some OMP directives.
if( DBS_CXX_IS_BULLSEYE OR ${CMAKE_CXX_COMPILER} MATCHES "[ps]pu-" )
   list( REMOVE_ITEM test_sources ${PROJECT_SOURCE_DIR}/tstOMP.cc )
endif()

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories( 
   ${PROJECT_SOURCE_DIR}      # headers for tests
   ${PROJECT_SOURCE_DIR}/..   # headers for package
   ${PROJECT_BINARY_DIR}/..   # config.h
   )

# ---------------------------------------------------------------------------- #
# Build library for test directory
# ---------------------------------------------------------------------------- #

add_component_library( 
   TARGET       Lib_c4_test
   TARGET_DEPS  Lib_c4
   LIBRARY_NAME c4_test 
   SOURCES      "${test_lib_sources}" )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #
set( test_deps
   Lib_c4
#   ${PAPI_LIBRARY}
#   ${omplib}
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

# This test also has special arguments and spawns mpirun indirectly so
# we register it as a scalar test.
set( special_test ${PROJECT_SOURCE_DIR}/tstApplicationUnitTest.cc )
add_scalar_tests(
   SOURCES   "${special_test}" 
   DEPS      "${test_deps}"
   TEST_ARGS "--np 1;--np 2;--np 4"
   APPLICATION_UNIT_TEST  )
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

if( DBS_CXX_IS_BULLSEYE OR ${CMAKE_CXX_COMPILER} MATCHES "[ps]pu-" )
   # no-op
else()
   # tstOMP should be run alone (not in paralle with other tests) because
   # it tries to use all available cores as OMP threads.
   # Also mark with 'nomemcheck' as valgrind can mess up the timing
   # checks. 
   if( ${DRACO_C4} MATCHES MPI )
      set_tests_properties( c4_tstOMP_1 
         PROPERTIES 
         RUN_SERIAL "ON"
      )
      set_tests_properties( c4_tstOMP_1 c4_tstOMP_2 c4_tstOMP_4 
         PROPERTIES 
         LABELS "nomemcheck"
      )
   endif()
  if( ${CMAKE_GENERATOR} MATCHES Xcode )
     SET_TARGET_PROPERTIES( Ut_c4_tstOMP_exe PROPERTIES XCODE_ATTRIBUTE_ENABLE_OPENMP_SUPPORT YES )
  endif()
endif()

 # RUN_SERIAL
 #       Do not run this test in parallel with any other test.
 #       Use this option in conjunction with the ctest_test PARALLEL_LEVEL
 #       option to specify that this test should not be run in parallel with
 #       any other tests.
