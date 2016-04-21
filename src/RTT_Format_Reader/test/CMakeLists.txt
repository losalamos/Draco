#-----------------------------*-cmake-*----------------------------------------#
# file   config/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for RTT_Format_Reader/test.
# note   Copyright (C) 2016, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#
project( RTT_Format_Reader_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB test_sources *.cc )
file( GLOB mesh_input_files *.mesh.in )

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories(
   ${PROJECT_SOURCE_DIR}      # headers for tests
   ${PROJECT_SOURCE_DIR}/..   # headers for package
   ${PROJECT_BINARY_DIR}/..   # config.h
)

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

set(test_deps Lib_RTT_Format_Reader Lib_mesh_element Lib_dsxx )
add_scalar_tests(
   SOURCES "${test_sources}"
   DEPS    "${test_deps}" )

# ---------------------------------------------------------------------------- #
# Copy input files to test directory (binary_dir)
# ---------------------------------------------------------------------------- #

provide_aux_files(
   FILES    "${mesh_input_files}"
   SRC_EXT  ".mesh.in"
   DEST_EXT ".mesh" )
