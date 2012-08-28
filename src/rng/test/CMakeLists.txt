# ---------------------------------------------------------------------------- #
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
# ---------------------------------------------------------------------------- #
# $Id$
# ---------------------------------------------------------------------------- #

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

add_component_library( Lib_rng_test ${PROJECT_NAME} ${test_lib_sources} )
add_dependencies( Lib_rng_test Lib_rng )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_scalar_tests(
   SOURCES "${test_sources}"
   DEPS    "Lib_dsxx;${GSL_LIBRARIES};" )

if( RANDOM123_FOUND)
   # Random123 unit tests almost all use assert to check for failures,
   # and CTest seems to identify the subsequent abort as a failure
   # regardless of the pass or fail regex.  ut_ars deviates from
   # this pattern by writing "error" plus an informative message to
   # stderr.
   add_scalar_tests(
      SOURCES    "${random123_unit_tests}"
      DEPS       "Lib_dsxx;${GSL_LIBRARIES};"
      PASS_REGEX "."
      FAIL_REGEX "failed;error" )

   add_scalar_tests(
      SOURCES    "${random123_known_answer_tests}"
      DEPS       "Lib_dsxx;"
      PASS_REGEX "PASSED"
      FAIL_REGEX "FAILED" )

   configure_file( ${PROJECT_SOURCE_DIR}/kat_vectors
      ${PROJECT_BINARY_DIR}/kat_vectors COPYONLY )
endif()
