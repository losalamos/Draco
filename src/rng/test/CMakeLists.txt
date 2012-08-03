# ---------------------------------------------------------------------------- #
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
# ---------------------------------------------------------------------------- #
# $Id$
# ---------------------------------------------------------------------------- #

project( rng_test C CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_lib_sources 
   ${PROJECT_SOURCE_DIR}/rng_test.cc 
   ${PROJECT_SOURCE_DIR}/rng_test.hh
)
file( GLOB test_sources *.cc )
list( REMOVE_ITEM test_sources ${test_lib_sources} )

set( random123_sources
   ${PROJECT_SOURCE_DIR}/tstTF_Gen.cc
   ${PROJECT_SOURCE_DIR}/tstTF_Rnd_Control_Inline.cc
)

if( RANDOM123_FOUND )
   file( GLOB random123_sources *.cpp *.c )
   list( APPEND test_sources ${random123_sources} )
else()
   list( REMOVE_ITEM test_sources ${random123_sources} )
endif()

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

add_component_library( Lib_rng_test ${PROJECT_NAME} ${test_lib_sources} )
add_dependencies( Lib_rng_test Lib_rng )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_scalar_tests(
   SOURCES "${test_sources}"
   DEPS    "Lib_dsxx;${GSL_LIBRARIES};" )

# ---------------------------------------------------------------------------- #
# Copy input files to test directory
# ---------------------------------------------------------------------------- #

if( RANDOM123_FOUND )
   configure_file( ${PROJECT_SOURCE_DIR}/kat_vectors
      ${PROJECT_BINARY_DIR}/kat_vectors COPYONLY )
endif()
