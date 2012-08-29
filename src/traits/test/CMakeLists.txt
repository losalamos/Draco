#-----------------------------*-cmake-*----------------------------------------#
# file   config/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for traits/test.
# note   Copyright (C) 2010-2012, Los Alamos National Security
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#
project( traits_test CXX )

file( GLOB test_sources *.cc )

# Find local headers
include_directories( 
   ${PROJECT_SOURCE_DIR}      # headers for tests
   ${PROJECT_SOURCE_DIR}/..   # headers for package
   ${dsxx_BINARY_DIR}         # headers for ds++
   ${draco_src_dir_SOURCE_DIR}
   ${MPI_INCLUDE_PATH} )

add_scalar_tests(
   SOURCES "${test_sources}"
   DEPS    "Lib_dsxx"
   )

