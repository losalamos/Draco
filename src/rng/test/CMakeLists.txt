#--------------------------------------------*-cmake-*---------------------------------------------#
# file   rng/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for rng/test.
# note   Copyright (C) 2016-2019, Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
project( rng_test C CXX )

#--------------------------------------------------------------------------------------------------#
# Source files
#--------------------------------------------------------------------------------------------------#

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
   ${PROJECT_SOURCE_DIR}/ut_M128.cpp
   ${PROJECT_SOURCE_DIR}/ut_uniform.cpp
)
# [2020-08-21 KT] The code in ut_gsl.c is not valid C++.  Because this test is a mix of C and C++,
# it causes issues. Most compilers are happy to compile this code with C, visual studio is not.
# When marked as C++, this test fails with:
#
# "C:\projects\draco\src\rng\test\ut_gsl.c(62,1): error C2086: 'const gsl_rng_type \
#    *gsl_rng_cbrng': redefinition"
#
# Indeed, reviewing the Random123 file gsl_cbrng.h reveals that that gsl_rng_type is declared twice,
# first on line 68 and then again on line 125 (with an assignment this time). We can continue
# building and running this tests on platforms that can processs it (all but MSVC right now).
if( NOT MSVC )
   list(APPEND random123_unit_tests ${PROJECT_SOURCE_DIR}/ut_gsl.c )
endif()
if( R123_USE_CXX11 )
   list( APPEND random123_unit_tests ${PROJECT_SOURCE_DIR}/ut_Engine.cpp )
endif()

# Random123 known-answer tests
set( random123_known_answer_tests
   ${PROJECT_SOURCE_DIR}/kat_c.c
   ${PROJECT_SOURCE_DIR}/kat_cpp.cpp
)

# Some compilers (esp. MSVC) have trouble with compiling these files as C11 code, so tell cmake
# to treat them as C++14. ut_gsl.c has a bug that prevents it from compiling as C++, so we omit it
# from this list.
if( MSVC )
  set_source_files_properties(
     ${PROJECT_SOURCE_DIR}/time_serial.c
     ${PROJECT_SOURCE_DIR}/ut_ars.c
     ${PROJECT_SOURCE_DIR}/kat_c.c
     PROPERTIES LANGUAGE CXX )
endif()

#--------------------------------------------------------------------------------------------------#
# Build Unit tests
#--------------------------------------------------------------------------------------------------#

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

 if( USE_CUDA AND NOT CMAKE_CXX_COMPILER_ID MATCHES XL )

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

#--------------------------------------------------------------------------------------------------#
# end rng/test/CMakeLists.txt
#--------------------------------------------------------------------------------------------------#
