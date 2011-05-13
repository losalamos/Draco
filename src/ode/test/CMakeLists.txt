# ---------------------------------------------------------------------------- #
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
# ---------------------------------------------------------------------------- #
# $Id$
# ---------------------------------------------------------------------------- #

project( ode_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_lib_sources 
   ${PROJECT_SOURCE_DIR}/ode_test.cc 
   ${PROJECT_SOURCE_DIR}/ode_test.hh
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

add_component_library( Lib_ode_test ${PROJECT_NAME} "${test_lib_sources}" )
add_dependencies( Lib_ode_test Lib_ode )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_scalar_tests(
   SOURCES "${test_sources}"
   DEPS    "Lib_c4;Lib_dsxx;${MPI_LIBRARY}"
   )

