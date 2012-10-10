#-----------------------------*-cmake-*----------------------------------------#
# file   config/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for rng/test.
# note   Copyright (C) 2010-2012, Los Alamos National Security
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#
project( rng_test C CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_lib_sources 
   ${PROJECT_SOURCE_DIR}/rng_test.cc 
   ${PROJECT_SOURCE_DIR}/rng_test.hh
)

set( test_sources
   ${PROJECT_SOURCE_DIR}/tstRnd_Control_Inline.cc
   ${PROJECT_SOURCE_DIR}/tstSubrandom_Sequence.cc
)

# Random123 unit tests
set( random123_unit_tests
   ${PROJECT_SOURCE_DIR}/ut_aes.cpp
   ${PROJECT_SOURCE_DIR}/ut_ars.c
   ${PROJECT_SOURCE_DIR}/ut_carray.cpp
   ${PROJECT_SOURCE_DIR}/ut_gsl.c
   ${PROJECT_SOURCE_DIR}/ut_Engine.cpp
   ${PROJECT_SOURCE_DIR}/ut_M128.cpp
)

# Random123 known-answer tests
set( random123_known_answer_tests
   ${PROJECT_SOURCE_DIR}/kat_c.c
   ${PROJECT_SOURCE_DIR}/kat_cpp.cpp
   ${PROJECT_SOURCE_DIR}/kat_u01_c.c
   ${PROJECT_SOURCE_DIR}/kat_u01_cpp.cpp
)

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

add_component_library( 
   TARGET       Lib_rng_test
   TARGET_DEPS  Lib_rng
   LIBRARY_NAME ${PROJECT_NAME} 
   SOURCES      "${test_lib_sources}" )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_scalar_tests( 
   SOURCES "${test_sources}"
   DEPS    "Lib_rng_test" ) # dsxx;${GSL_LIBRARIES};" )

if( RANDOM123_FOUND)
   # Random123 unit tests almost all use assert to check for failures,
   # and CTest seems to identify the subsequent abort as a failure
   # regardless of the pass or fail regex.  ut_ars deviates from
   # this pattern by writing "error" plus an informative message to
   # stderr.
   #
   # Because assert is removed at compile time when NDEBUG is defined,
   # these tests should be run only with Debug and RelWithDebInfo
   # builds.
   string( TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UPPER )
   if( "${CMAKE_BUILD_TYPE_UPPER}" STREQUAL "DEBUG" OR
         "${CMAKE_BUILD_TYPE_UPPER}" STREQUAL "RELWITHDEBINFO" )
      add_scalar_tests(
         SOURCES    "${random123_unit_tests}"
         PASS_REGEX "."
         FAIL_REGEX "failed;error;Failure" )
   endif()

   add_scalar_tests(
      SOURCES    "${random123_known_answer_tests}"
      PASS_REGEX "PASSED"
      FAIL_REGEX "FAILED" )

   if( USE_CUDA AND EXISTS "${CUDA_NVCC_EXECUTABLE}" )
      cuda_add_executable( kat_cuda
         ${PROJECT_SOURCE_DIR}/kat_cuda.cu )

      add_test(
         NAME    rng_kat_cuda
         COMMAND $<TARGET_FILE:kat_cuda> )

      set_tests_properties( rng_kat_cuda
         PROPERTIES
         PASS_REGULAR_EXPRESSION "PASSED"
         FAIL_REGULAR_EXPRESSION "FAILED" )

      cuda_add_executable( kat_u01_cuda
         ${PROJECT_SOURCE_DIR}/kat_u01_cuda.cu )

      add_test(
         NAME    rng_kat_u01_cuda
         COMMAND $<TARGET_FILE:kat_u01_cuda> )

      set_tests_properties( rng_kat_u01_cuda
         PROPERTIES
         PASS_REGULAR_EXPRESSION "PASSED"
         FAIL_REGULAR_EXPRESSION "FAILED" )
   endif()

   configure_file( ${PROJECT_SOURCE_DIR}/kat_vectors
      ${PROJECT_BINARY_DIR}/kat_vectors COPYONLY )
endif()
