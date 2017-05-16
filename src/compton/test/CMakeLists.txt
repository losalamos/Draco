#-----------------------------*-cmake-*----------------------------------------#
# file   compton/test/CMakeLists.txt
# author Kendra Keady
# date   2017 February 28
# brief  Generate build project files for compton/test.
# note   Copyright (C) 2017, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

project( compton_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

 set( test_sources
  ${PROJECT_SOURCE_DIR}/tCompton.cc )
 file( GLOB compton_data *.compton )

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories(
  ${PROJECT_SOURCE_DIR}      # headers for tests
  ${PROJECT_BINARY_DIR}      # test/config.h
  ${PROJECT_SOURCE_DIR}/..   # headers for package
  )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_scalar_tests(
  SOURCES "${test_sources}"
  DEPS    "Lib_compton" )

#------------------------------------------------------------------------------#
# Copy data files
#------------------------------------------------------------------------------#

provide_aux_files( FILES "${compton_data}" )

# ---------------------------------------------------------------------------- #
# end of CMakeLists.txt
# ---------------------------------------------------------------------------- #
