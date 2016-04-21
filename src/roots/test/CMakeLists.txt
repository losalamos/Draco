#-----------------------------*-cmake-*----------------------------------------#
# file   roots/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for roots/test.
# note   Copyright (C) 2016, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#
project( roots_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB test_sources *.cc )
file( GLOB mesh_input_files *.in )
file( GLOB python_scripts python/*.py )
file( GLOB bench_files bench/*.0001 )

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories(
   ${PROJECT_SOURCE_DIR}      # headers for tests
   ${PROJECT_SOURCE_DIR}/..   # headers for package
   ${c4_BINARY_DIR}           # c4/config.h
   ${MPI_CXX_INCLUDE_PATH}
   )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_scalar_tests(
   SOURCES "${test_sources}"
   DEPS "Lib_roots;Lib_linear" )
