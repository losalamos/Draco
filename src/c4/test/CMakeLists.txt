#-----------------------------*-cmake-*----------------------------------------#
# file   c4/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for c4/test.
# note   Copyright (C) 2016-2019, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
project( c4_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB test_sources *.cc )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

# tstAbort doesn't really work as a unit test because it always calls MPI_abort.
# This causes an uncatchable error condition.  I suppose we could wrap the test
# in a python or cmake script - maybe we can do this sometime in the future.
list( REMOVE_ITEM test_sources ${PROJECT_SOURCE_DIR}/tstAbort.cc )

# This test has special arguments and special pass/fail criteria.
set( special_test ${PROJECT_SOURCE_DIR}/phw.cc )
add_parallel_tests(
   SOURCES   "${special_test}"
   DEPS      Lib_c4
   PE_LIST   "1;2;4"
   TEST_ARGS "hello"
   PASS_REGEX "Hello"
   FAIL_REGEX "cruel"
   )
list( REMOVE_ITEM test_sources ${special_test} )

# This tests runs with extra arguments
set( special_test ${PROJECT_SOURCE_DIR}/tstParallelUnitTest.cc )
add_parallel_tests(
   SOURCES   "${special_test}"
   DEPS      Lib_c4
   PE_LIST   "1;2;4"
   TEST_ARGS "--version;--dummy; "
   PASS_REGEX "tstParallelUnitTest[.ex]*: version"
   )
list( REMOVE_ITEM test_sources ${special_test} )

# This test uses MPI+OMP:
# - Bullseye doesn't like some OMP directives.
# - Timing tests are broken if the test has to be pushed to a virtual
#   machine (aprun).
set( special_test ${PROJECT_SOURCE_DIR}/tstOMP.cc )
list( REMOVE_ITEM test_sources ${special_test} )
if( NOT OPENMP_FOUND )
   # no-op
elseif( DBS_CXX_IS_BULLSEYE )
   # no-op
else()
   add_parallel_tests(
      SOURCES   "${special_test}"
      DEPS      Lib_c4
      PE_LIST   "1;2"
      MPI_PLUS_OMP
      )
endif()

# Tests that expect exactly 2 PE (omit for DRACO_C4=SCALAR builds):
set( two_pe_tests ${PROJECT_SOURCE_DIR}/tstPingPong.cc )
list( APPEND two_pe_tests ${PROJECT_SOURCE_DIR}/tstC4_Status.cc )
list( REMOVE_ITEM test_sources ${two_pe_tests} )
if( NOT ${DRACO_C4} STREQUAL "SCALAR" )
  add_parallel_tests(
    SOURCES   "${two_pe_tests}"
    DEPS      Lib_c4
    PE_LIST   "2"
    )
endif()

# Exclude tstTime from valgrind checks.
set( nomemcheck_tests
  ${PROJECT_SOURCE_DIR}/tstTime.cc
  ${PROJECT_SOURCE_DIR}/tstGlobalTimer.cc
  )
list( REMOVE_ITEM test_sources ${nomemcheck_tests} )
add_parallel_tests(
  SOURCES   "${nomemcheck_tests}"
  DEPS      Lib_c4
  PE_LIST   "1;2;4"
  LABEL     "nomemcheck"
  )

# The remainder of the tests
add_parallel_tests(
   SOURCES   "${test_sources}"
   DEPS      Lib_c4
   PE_LIST   "1;2;4"
   )

#------------------------------------------------------------------------------#
# Python based test runners
#------------------------------------------------------------------------------#
if( NOT APPLE )
  include( ApplicationUnitTest )
  add_app_unit_test(
    DRIVER    ${CMAKE_CURRENT_SOURCE_DIR}/tstXthi.py
    APP       $<TARGET_FILE_DIR:Exe_xthi>/$<TARGET_FILE_NAME:Exe_xthi>
    PE_LIST   "1;2" )
  add_app_unit_test(
    DRIVER    ${CMAKE_CURRENT_SOURCE_DIR}/tstYthi.py
    APP       $<TARGET_FILE_DIR:Exe_ythi>/$<TARGET_FILE_NAME:Exe_ythi>
    PE_LIST   "1;2" )
endif()

#------------------------------------------------------------------------------#
# End src/c4/test/CMakeLists.txt
#------------------------------------------------------------------------------#
