#-----------------------------*-cmake-*----------------------------------------#
# file   config/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for rng/test.
# note   Copyright (C) 2016-2017, Los Alamos National Security, LLC.
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

if( USE_CUDA AND EXISTS "${CUDA_NVCC_EXECUTABLE}" )
  cuda_add_executable( kat_cuda
    ${PROJECT_SOURCE_DIR}/kat_cuda.cu )

  add_test(
    NAME    rng_kat_cuda
    COMMAND $<TARGET_FILE:kat_cuda> )

  set_tests_properties( rng_kat_cuda PROPERTIES
    PASS_REGULAR_EXPRESSION "PASSED"
    FAIL_REGULAR_EXPRESSION "FAILED" )
endif()

add_component_executable(
  TARGET  Exe_time_serial
  SOURCES ${PROJECT_SOURCE_DIR}/time_serial.c
  NOEXPORT
  )

configure_file( ${PROJECT_SOURCE_DIR}/kat_vectors
  ${PROJECT_BINARY_DIR}/kat_vectors COPYONLY )

# ---------------------------------------------------------------------------- #
# end rng/test/CMakeLists.txt
# ---------------------------------------------------------------------------- #
