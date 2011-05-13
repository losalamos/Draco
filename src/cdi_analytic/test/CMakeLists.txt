# ---------------------------------------------------------------------------- #
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
# ---------------------------------------------------------------------------- #
# $Id$
# ---------------------------------------------------------------------------- #

project( cdi_analytic_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_lib_sources 
   ${PROJECT_SOURCE_DIR}/cdi_analytic_test.cc 
   ${PROJECT_SOURCE_DIR}/cdi_analytic_test.hh
)
file( GLOB test_sources *.cc )
list( REMOVE_ITEM test_sources ${test_lib_sources} )

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories( 
   ${PROJECT_SOURCE_DIR}      # headers for tests
   ${PROJECT_SOURCE_DIR}/..   # headers for package
   ${PROJECT_BINARY_DIR}/..   # c4/config.h
   ${MPI_INCLUDE_PATH} 
)

# ---------------------------------------------------------------------------- #
# Build library for test directory
# ---------------------------------------------------------------------------- #

add_component_library( Lib_cdi_analytic_test ${PROJECT_NAME} ${test_lib_sources} )
add_dependencies(      Lib_cdi_analytic_test Lib_cdi_analytic )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_scalar_tests(
   SOURCES "${test_sources}" 
   DEPS    "Lib_parser;Lib_cdi;Lib_c4;${MPI_LIBRARIES};Lib_dsxx" )




