# ---------------------------------------------------------------------------- #
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
# ---------------------------------------------------------------------------- #
# $Id$
# ---------------------------------------------------------------------------- #

project( dsxx_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_lib_sources 
   ${PROJECT_SOURCE_DIR}/ds_test.cc 
   ${PROJECT_SOURCE_DIR}/ds_test.hh
)
file( GLOB test_sources *.cc )
list( REMOVE_ITEM test_sources ${test_lib_sources} )

# BUG: artf21615 (tf.lanl.gov): Safe_Ptr has a design defect
list( REMOVE_ITEM test_sources 
   ${PROJECT_SOURCE_DIR}/tstSafe_Ptr.cc
   ${PROJECT_SOURCE_DIR}/tstDBC_Ptr.cc )

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories( 
   ${PROJECT_SOURCE_DIR}      # headers for tests
   ${PROJECT_SOURCE_DIR}/..   # headers for package
   ${PROJECT_BINARY_DIR}/.. ) # config.h

# ---------------------------------------------------------------------------- #
# Build library for test directory
# ---------------------------------------------------------------------------- #

add_component_library( Lib_dsxx_test ds++_test     "${test_lib_sources}" )
#                      ^target_name  ^library_name ^source list
add_dependencies( Lib_dsxx_test Lib_dsxx )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_scalar_tests( 
   SOURCES "${test_sources}" 
   )



