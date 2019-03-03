#-----------------------------*-cmake-*----------------------------------------#
# file   timestep/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for timestep/test.
# note   Copyright (C) 2016-2019, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
project( timestep_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_lib_sources
   ${PROJECT_SOURCE_DIR}/dummy_package.cc
   ${PROJECT_SOURCE_DIR}/dummy_package.hh
)
file( GLOB test_sources *.cc )
list( REMOVE_ITEM test_sources ${test_lib_sources} )

# ---------------------------------------------------------------------------- #
# Build library for test directory
# ---------------------------------------------------------------------------- #

add_component_library(
   TARGET       Lib_timestep_test
   TARGET_DEPS  Lib_timestep
   LIBRARY_NAME ${PROJECT_NAME}
   SOURCES      "${test_lib_sources}" )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_parallel_tests(
   SOURCES "${test_sources}"
   PE_LIST "1;2;4"
   DEPS    "Lib_timestep_test;${MPI_CXX_LIBRARIES};${PAPI_LIBRARY}" )

#------------------------------------------------------------------------------#
# End timestep/test/CMakeLists.txt
#------------------------------------------------------------------------------#
