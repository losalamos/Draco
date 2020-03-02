#-----------------------------*-cmake-*----------------------------------------#
# file   rng/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for rng/test.
# note   Copyright (C) 2016-2019, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
project( rng_test C CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_sources
   ${PROJECT_SOURCE_DIR}/tstRnd_Control_Inline.cc
   ${PROJECT_SOURCE_DIR}/tstSubrandom_Sequence.cc
   ${PROJECT_SOURCE_DIR}/tstCounter_RNG.cc
)

# Random123 unit tests (these tests have special PASS/FAIL REGEX conditions)
set( random123_unit_tests
   ${PROJECT_SOURCE_DIR}/ut_aes.cpp
   ${PROJECT_SOURCE_DIR}/ut_ars.c
   ${PROJECT_SOURCE_DIR}/ut_carray.cpp
   ${PROJECT_SOURCE_DIR}/ut_gsl.c
   ${PROJECT_SOURCE_DIR}/ut_M128.cpp
   ${PROJECT_SOURCE_DIR}/ut_uniform.cpp
)

if( R123_USE_CXX11 )
   list( APPEND random123_unit_tests
      ${PROJECT_SOURCE_DIR}/ut_Engine.cpp )
endif()

# Random123 known-answer tests
set( random123_known_answer_tests
   ${PROJECT_SOURCE_DIR}/kat_c.c
   ${PROJECT_SOURCE_DIR}/kat_cpp.cpp
)

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_scalar_tests(
   SOURCES "${test_sources}"
   DEPS    "Lib_rng" )

add_scalar_tests(
   SOURCES    "${random123_unit_tests}"
   DEPS       "Lib_rng"
   PASS_REGEX "OK|SUCCESS"
   FAIL_REGEX "failed;error;Failure" )

add_scalar_tests(
  SOURCES    "${random123_known_answer_tests}"
  DEPS       "Lib_rng"
  PASS_REGEX "PASSED"
  FAIL_REGEX "FAILED" )
target_include_directories( Ut_rng_ut_uniform_exe
  PRIVATE $<BUILD_INTERFACE:${rng_SOURCE_DIR}> )

 if( HAVE_CUDA )

  set( test_sources_cuda
     ${PROJECT_SOURCE_DIR}/kat_cuda.cu)

  set( test_deps_cuda "Lib_rng" )

   add_scalar_tests(
      SOURCES  "${test_sources_cuda}"
      DEPS     "${test_deps_cuda}"
      PASS_REGEX "PASSED"
      FAIL_REGEX "FAILED" )

 endif()

add_component_executable(
  TARGET  Exe_time_serial
  SOURCES ${PROJECT_SOURCE_DIR}/time_serial.c
  NOEXPORT
  )
target_include_directories( Exe_time_serial
  PRIVATE
    $<BUILD_INTERFACE:${rng_BINARY_DIR}>
    $<BUILD_INTERFACE:${dsxx_BINARY_DIR}>
    ${RANDOM123_INCLUDE_DIR} )

configure_file( ${PROJECT_SOURCE_DIR}/kat_vectors
  ${PROJECT_BINARY_DIR}/kat_vectors COPYONLY )

# ---------------------------------------------------------------------------- #
# end rng/test/CMakeLists.txt
# ---------------------------------------------------------------------------- #
