#=============================================================================
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
#=============================================================================
# $Id$
#=============================================================================

project( pcgWrap_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_lib_sources 
	${PROJECT_SOURCE_DIR}/test_utils.cc 
	${PROJECT_SOURCE_DIR}/test_utils.hh)
file( GLOB test_sources *.cc )
list( REMOVE_ITEM test_sources ${test_lib_sources} )

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #
include_directories( 
   ${PROJECT_SOURCE_DIR}      # headers for tests
   ${PROJECT_SOURCE_DIR}/..   # headers for package
   ${PROJECT_BINARY_DIR}/..   # config.h
   )

# ---------------------------------------------------------------------------- #
# Build library for test directory
# ---------------------------------------------------------------------------- #

add_component_library( Lib_pcgWrap_test ${PROJECT_NAME} "${test_lib_sources}" )
add_dependencies( Lib_pcgWrap_test Lib_pcgWrap )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

set( test_deps 
   Lib_dsxx 
   ${PCG_LIBRARY} 
   )

add_scalar_tests(
   SOURCES "${test_sources}"
   DEPS    "${test_deps}" ) 

# ---------------------------------------------------------------------------- #
# Copy input files to test directory (binary_dir)
# ---------------------------------------------------------------------------- #

if( WIN32 )
   provide_aux_files( FILES "${GSL_DLL_LIBRARIES}" )
endif()

