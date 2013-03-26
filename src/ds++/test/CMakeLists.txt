#-----------------------------*-cmake-*----------------------------------------#
# file   src/ds++/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 April 28
# brief  Instructions for building ds++/test level Makefile.
# note   Copyright (C) 2010-2013 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#
project( dsxx_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB test_sources *.cc )

# BUG: artf21615 (tf.lanl.gov): Safe_Ptr has a design defect
list( REMOVE_ITEM test_sources 
   ${PROJECT_SOURCE_DIR}/tstSafe_Ptr.cc
   ${PROJECT_SOURCE_DIR}/tstDBC_Ptr.cc )

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories( 
   ${PROJECT_SOURCE_DIR}      # headers for tests
   ${PROJECT_SOURCE_DIR}/..   # headers for package
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
