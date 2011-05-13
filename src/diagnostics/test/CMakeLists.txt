# ---------------------------------------------------------------------------- #
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
# ---------------------------------------------------------------------------- #
# $Id$
# ---------------------------------------------------------------------------- #

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
   )

# ---------------------------------------------------------------------------- #
# Build library for test directory
# ---------------------------------------------------------------------------- #

add_component_library( Lib_diagnostics_test ${PROJECT_NAME} "${test_lib_sources}" )
add_dependencies(      Lib_diagnostics_test Lib_diagnostics Lib_c4 )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_parallel_tests(
   SOURCES "${test_sources}"
   PE_LIST "1"
   DEPS    "Lib_c4;Lib_dsxx;${MPI_LIBRARIES}" )




