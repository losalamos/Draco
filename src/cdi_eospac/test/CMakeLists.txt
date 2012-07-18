# ---------------------------------------------------------------------------- #
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
# ---------------------------------------------------------------------------- #
# $Id: CMakeLists.txt 6292 2011-12-07 01:40:03Z kellyt $
# ---------------------------------------------------------------------------- #

project( cdi_eospac_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_lib_sources 
   ${PROJECT_SOURCE_DIR}/cdi_eospac_test.cc 
   ${PROJECT_SOURCE_DIR}/cdi_eospac_test.hh
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

add_component_library( Lib_cdi_eospac_test ${PROJECT_NAME} ${test_lib_sources} )
add_dependencies(      Lib_cdi_eospac_test Lib_cdi_eospac )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_scalar_tests(
   SOURCES "${test_sources}" 
   DEPS    "Lib_cdi;Lib_c4;${MPI_LIBRARIES};Lib_dsxx;${PAPI_LIBRARY}" )




