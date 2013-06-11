#-----------------------------*-cmake-*----------------------------------------#
# file   config/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for meshReaders/test.
# note   Copyright (C) 2010-2013, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#
project( meshReaders_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_lib_sources 
   ${PROJECT_SOURCE_DIR}/meshReaders_test.cc 
   ${PROJECT_SOURCE_DIR}/meshReaders_test.hh
)
file( GLOB test_sources *.cc )
file( GLOB mesh_input_files *.mesh.in )
list( REMOVE_ITEM test_sources ${test_lib_sources} )

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
   TARGET       Lib_meshReaders_test
   TARGET_DEPS  Lib_meshReaders
   LIBRARY_NAME ${PROJECT_NAME} 
   SOURCES "${test_lib_sources}" )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_scalar_tests(
   SOURCES "${test_sources}"
   DEPS    "Lib_dsxx" )

# ---------------------------------------------------------------------------- #
# Copy mesh files to test directory (binary_dir)
# ---------------------------------------------------------------------------- #

provide_aux_files(
   FILES    "${mesh_input_files}"
   SRC_EXT  ".mesh.in"
   DEST_EXT ".mesh" )


