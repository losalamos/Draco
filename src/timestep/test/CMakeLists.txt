# ---------------------------------------------------------------------------- #
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
# ---------------------------------------------------------------------------- #
# $Id$
# ---------------------------------------------------------------------------- #

project( timestep_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_lib_sources 
   ${PROJECT_SOURCE_DIR}/timestep_test.cc 
   ${PROJECT_SOURCE_DIR}/timestep_test.hh
   ${PROJECT_SOURCE_DIR}/dummy_package.cc 
   ${PROJECT_SOURCE_DIR}/dummy_package.hh
)
file( GLOB test_sources *.cc )
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

add_component_library( Lib_timestep_test ${PROJECT_NAME} "${test_lib_sources}" )
add_dependencies( Lib_timestep_test Lib_timestep )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_parallel_tests(
   SOURCES "${test_sources}"
   PE_LIST "1;2;4"
   DEPS    "Lib_c4;Lib_dsxx;${MPI_LIBRARIES};${INTEL_LIBM}" )

