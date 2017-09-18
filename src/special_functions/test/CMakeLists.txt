#-----------------------------*-cmake-*----------------------------------------#
# file   special_functions/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 Aug 1
# brief  Generate build project files for special_functions/test
# note   Copyright (C) 2016, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.9.0)
project( special_functions_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #
file( GLOB test_sources *.cc )

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories(
   ${PROJECT_SOURCE_DIR}      # headers for tests
   ${PROJECT_SOURCE_DIR}/..   # headers for package
   ${ode_BINARY_DIR}          # ode/config.h
   ${diagnostics_BINARY_DIR}  # diagnostics/config.h
   )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_scalar_tests(
   SOURCES "${test_sources}"
   DEPS    "Lib_special_functions;Lib_diagnostics" )

# ---------------------------------------------------------------------------- #
# Copy vendor dll files to test directory (binary_dir)
# ---------------------------------------------------------------------------- #

if( WIN32 )
   provide_aux_files( FILES "${GSL_DLL_LIBRARIES}" )
endif()
