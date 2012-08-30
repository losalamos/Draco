#-----------------------------*-cmake-*----------------------------------------#
# file   config/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for diagnostics/test.
# note   Copyright (C) 2010-2012, Los Alamos National Security
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#
project( diagnostics_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_lib_sources 
   ${PROJECT_SOURCE_DIR}/diagnostics_test.cc 
   ${PROJECT_SOURCE_DIR}/diagnostics_test.hh
)
file( GLOB test_sources *.cc )
list( REMOVE_ITEM test_sources ${test_lib_sources} )

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories( 
   ${PROJECT_SOURCE_DIR}      # headers for tests
   ${PROJECT_SOURCE_DIR}/..   # headers for package
   ${PROJECT_BINARY_DIR}/..   # config.h
   ${c4_BINARY_DIR}           # c4/config.h
   ${MPI_CXX_INCLUDE_PATH}
   )

# ---------------------------------------------------------------------------- #
# Build library for test directory
# ---------------------------------------------------------------------------- #

add_component_library( 
   TARGET       Lib_diagnostics_test
   TARGET_DEPS  Lib_diagnostics
   LIBRARY_NAME ${PROJECT_NAME}
   SOURCES      "${test_lib_sources}" )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_parallel_tests(
   SOURCES "${test_sources}"
   PE_LIST "1"
   DEPS    "Lib_c4;Lib_dsxx;${MPI_LIBRARIES};${PAPI_LIBRARY}" )




