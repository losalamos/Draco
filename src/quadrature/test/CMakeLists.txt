#-----------------------------*-cmake-*----------------------------------------#
# file   quadrature/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for quadrature/test.
# note   Copyright (C) 2016-2019, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
project( quadrature_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_lib_sources
  ${PROJECT_SOURCE_DIR}/quadrature_test.cc
  ${PROJECT_SOURCE_DIR}/quadrature_test.hh
  )
file( GLOB test_sources *.cc )
list( REMOVE_ITEM test_sources ${test_lib_sources} )
file( GLOB input_files *.inp )

# ---------------------------------------------------------------------------- #
# Build library for test directory
# ---------------------------------------------------------------------------- #

add_component_library(
   TARGET       Lib_quadrature_test
   TARGET_DEPS  Lib_quadrature
   LIBRARY_NAME quadrature_test
   SOURCES      "${test_lib_sources}"
   PREFIX       Draco )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_scalar_tests(
  SOURCES "${test_sources}"
  DEPS    Lib_quadrature_test )

# ---------------------------------------------------------------------------- #
# Copy input files to test directory (binary_dir)
# ---------------------------------------------------------------------------- #

if( WIN32 )
   provide_aux_files( FILES "${GSL_DLL_LIBRARIES}" )
endif()

#------------------------------------------------------------------------------#
# End quadrature/test/CMakeLists.txt
#------------------------------------------------------------------------------#
