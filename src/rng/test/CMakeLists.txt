#-----------------------------*-cmake-*----------------------------------------#
# file   config/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for rng/test.
# note   Copyright (C) 2010-2013, Los Alamos National Security
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#
project( rng_test C CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_sources
   ${PROJECT_SOURCE_DIR}/tstRnd_Control_Inline.cc
   ${PROJECT_SOURCE_DIR}/tstSubrandom_Sequence.cc
)

if( RANDOM123_FOUND )
   list( APPEND test_sources ${PROJECT_SOURCE_DIR}/tstCounter_RNG.cc )
endif()

# Random123 unit tests
set( random123_unit_tests
   ${PROJECT_SOURCE_DIR}/ut_aes.cpp
   ${PROJECT_SOURCE_DIR}/ut_ars.c
   ${PROJECT_SOURCE_DIR}/ut_carray.cpp
   ${PROJECT_SOURCE_DIR}/ut_gsl.c
   ${PROJECT_SOURCE_DIR}/ut_M128.cpp
   ${PROJECT_SOURCE_DIR}/ut_uniform.cpp
)

if( R123_USE_CXX11 ) #if( NOT WIN32 )
   # Visual Studio 2010 does not support uniform_int_distribution<int> from C++11/C++0x.
   list( APPEND random123_unit_tests 
      ${PROJECT_SOURCE_DIR}/ut_Engine.cpp )
endif()

# Random123 known-answer tests
set( random123_known_answer_tests
   ${PROJECT_SOURCE_DIR}/kat_c.c
   ${PROJECT_SOURCE_DIR}/kat_cpp.cpp
)

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories( 
   ${PROJECT_SOURCE_DIR}      # headers for tests
   ${PROJECT_SOURCE_DIR}/..   # headers for package
)

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_scalar_tests( 
   SOURCES "${test_sources}"
   DEPS    "Lib_rng;${GSL_LIBRARIES}" )
   
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
         DEPS       "Lib_rng"
         PASS_REGEX "OK|SUCCESS"
         FAIL_REGEX "failed;error;Failure" )
   endif()

   add_scalar_tests(
      SOURCES    "${random123_known_answer_tests}"
      DEPS       "Lib_rng"
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
   endif()

   add_executable(time_serial
      ${PROJECT_SOURCE_DIR}/time_serial.c )
   set_target_properties( time_serial
     PROPERTIES 
       #OUTPUT_NAME time_serial
       PROJECT_LABEL time_serial
       FOLDER rng_test )

   configure_file( ${PROJECT_SOURCE_DIR}/kat_vectors
      ${PROJECT_BINARY_DIR}/kat_vectors COPYONLY )
endif()
