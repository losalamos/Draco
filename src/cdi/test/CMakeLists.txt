# ---------------------------------------------------------------------------- #
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
# ---------------------------------------------------------------------------- #
# $Id$
# ---------------------------------------------------------------------------- #

project( cdi_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_lib_sources 
   ${PROJECT_SOURCE_DIR}/cdi_test.cc 
   ${PROJECT_SOURCE_DIR}/cdi_test.hh
   ${PROJECT_SOURCE_DIR}/DummyEoS.cc
   ${PROJECT_SOURCE_DIR}/DummyGrayOpacity.cc
   ${PROJECT_SOURCE_DIR}/DummyMultigroupOpacity.cc
   ${PROJECT_SOURCE_DIR}/DummyOdfmgOpacity.cc   
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

add_component_library( Lib_cdi_test cdi_test "${test_lib_sources}" )
add_dependencies( Lib_cdi_test Lib_cdi )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_scalar_tests(
   SOURCES   "${test_sources}" 
   DEPS      "Lib_dsxx;${MPI_LIBRARIES};${INTEL_LIBM}" )

# source_group( cdi\\test FILES ${test_lib_sources} ${test_sources} )


