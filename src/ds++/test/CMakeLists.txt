#-----------------------------*-cmake-*----------------------------------------#
# file   src/ds++/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 April 28
# brief  Instructions for building ds++/test level Makefile.
# note   Copyright (C) 2010-2012 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#
project( dsxx_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_lib_sources 
   ${PROJECT_SOURCE_DIR}/ds_test.cc 
   ${PROJECT_SOURCE_DIR}/ds_test.hh
)
file( GLOB test_sources *.cc )
list( REMOVE_ITEM test_sources ${test_lib_sources} )

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
# Build library for test directory
# ---------------------------------------------------------------------------- #

add_component_library( 
   TARGET       Lib_dsxx_test
   TARGET_DEPS  Lib_dsxx
   LIBRARY_NAME ds++_test 
   SOURCES      "${test_lib_sources}" )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_scalar_tests( 
   SOURCES "${test_sources}" 
   )



