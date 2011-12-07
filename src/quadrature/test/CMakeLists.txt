#=============================================================================
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
#=============================================================================
# $Id$
#=============================================================================

project( quadrature_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_lib_sources 
	${PROJECT_SOURCE_DIR}/quadrature_test.cc 
	${PROJECT_SOURCE_DIR}/quadrature_test.hh)
file( GLOB test_sources *.cc )
list( REMOVE_ITEM test_sources ${test_lib_sources} )
file( GLOB input_files *.inp )

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #
include_directories( 
   ${PROJECT_SOURCE_DIR}      # headers for tests
   ${PROJECT_SOURCE_DIR}/..   # headers for package
   ${PROJECT_BINARY_DIR}/..   # config.h
   ${MPI_INCLUDE_PATH} )

# ---------------------------------------------------------------------------- #
# Build library for test directory
# ---------------------------------------------------------------------------- #

add_component_library( Lib_quadrature_test ${PROJECT_NAME} "${test_lib_sources}" )
add_dependencies( Lib_quadrature_test Lib_quadrature )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

set( test_deps
   Lib_parser
   Lib_c4
   Lib_special_functions
   Lib_ode
   Lib_dsxx
   ${GSL_LIBRARIES}
   ${MPI_LIBRARIES} 
   ${PAPI_LIBRARY} )

add_scalar_tests(
   SOURCES "${test_sources}"
   DEPS    "${test_deps}" ) 

# ---------------------------------------------------------------------------- #
# Copy input files to test directory (binary_dir)
# ---------------------------------------------------------------------------- #

if( WIN32 )
   provide_aux_files( FILES "${GSL_DLL_LIBRARIES}" )
endif()

