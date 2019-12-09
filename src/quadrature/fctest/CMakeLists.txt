#-----------------------------*-cmake-*----------------------------------------#
# file   quadrature/fctest/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2016 May 24
# brief  Generate build project files for quadrature/fctest
# note   Copyright (C) 2016, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
project( quadrature_fctest CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB test_sources *.cc )
#list( REMOVE_ITEM test_sources ${test_lib_sources} )
#file( GLOB input_files *.inp )

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

# Add the Fortran test and link against the f90 library from ftest

add_scalar_tests(
  SOURCES "${PROJECT_SOURCE_DIR}/ftstquadrature_interfaces.cc"
  DEPS    "Lib_quadrature_ftest;Lib_quadrature_test" )

# ---------------------------------------------------------------------------- #
# end quadrature/fctest/CMakeLists.txt
# ---------------------------------------------------------------------------- #

