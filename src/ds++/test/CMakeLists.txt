#-----------------------------*-cmake-*----------------------------------------#
# file   src/ds++/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 April 28
# brief  Instructions for building ds++/test level Makefile.
# note   Copyright (C) 2010-2015 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#
project( dsxx_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB test_sources *.cc )

list( REMOVE_ITEM test_sources
   ${PROJECT_SOURCE_DIR}/do_exception.cc )

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories(
   ${PROJECT_SOURCE_DIR}      # headers for tests
   ${Draco_SOURCE_DIR}/src    # headers for package
   ${PROJECT_BINARY_DIR}/.. ) # config.h

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

# Run ScalarUnitTest in 2 modes: (a) without additional command line
# arguments and (b) a special test to demonstrate the optional
# '--version' command line argument
add_scalar_tests(
   SOURCES ${PROJECT_SOURCE_DIR}/tstScalarUnitTest.cc
   DEPS    Lib_dsxx
   # Run first case without arguments
   # Run second case with the '--version' argument.
   TEST_ARGS " ;--version"
   #         ^^^ this space before the semicolon is important!
   # The 2nd case does not print "PASS" so we register the following
   # regex as PASS message.
   PASS_REGEX "tstScalarUnitTest[.ex]*: version"
)
list( REMOVE_ITEM test_sources
   ${PROJECT_SOURCE_DIR}/tstScalarUnitTest.cc )

# The remaining tests
add_scalar_tests(
   SOURCES "${test_sources}"
   DEPS    Lib_dsxx
   )

# (2013-11-13 KT):: tstdbc intentionally divides by zero (to test the
# isinfinity function). PGI issues a warning for this case that we
# don't care about.  Since PGI does not appear to use pragmas to
# supress warnings, I am turning off warnings for Ut_dsxx_tstdbc_exe
# when the compiler is PGI.
if( draco_isPGI )
   set_source_files_properties(
      ${PROJECT_SOURCE_DIR}/tstdbc.cc
      PROPERTIES
        COMPILE_FLAGS --no_warnings
   )
endif()

# ---------------------------------------------------------------------------- #
# This tests requries python.
# Python is not available on Darwin's MIC backends.
# ---------------------------------------------------------------------------- #
if( PYTHONINTERP_FOUND AND NOT HAVE_MIC )

  add_component_executable(
    TARGET      Exe_do_exception
    TARGET_DEPS Lib_dsxx
    SOURCES     ${PROJECT_SOURCE_DIR}/do_exception.cc
    NOEXPORT
    )
  add_test(
    NAME    dsxx_fpe_trap_test
    COMMAND ${PYTHON_EXECUTABLE}
            ${PROJECT_SOURCE_DIR}/python/test_fpe_trap.py
            $<TARGET_FILE:Exe_do_exception>
    )
  set_tests_properties( dsxx_fpe_trap_test
    PROPERTIES
      PASS_REGULAR_EXPRESSION ".*PASSED"
      FAIL_REGULAR_EXPRESSION ".*FAILED"
      LABELS                  "nomemcheck"
    )

endif() # PYTHONINTERP_FOUND
