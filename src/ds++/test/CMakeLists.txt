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

foreach( file ${test_sources} )
  get_filename_component( testname ${file} NAME_WE )
  add_executable( Ut_dsxx_${testname}_exe ${file})
  set_target_properties( 
     Ut_dsxx_${testname}_exe 
     PROPERTIES 
       OUTPUT_NAME ${testname} 
       VS_KEYWORD ${testname}
       PROJECT_LABEL Ut_dsxx
  )
  target_link_libraries( Ut_dsxx_${testname}_exe Lib_dsxx Lib_dsxx_test )
  add_test( dsxx_${testname} ${testname} )
  set_tests_properties( dsxx_${testname} PROPERTIES	
     PASS_REGULAR_EXPRESSION ".*[Tt]est: PASSED" 
     FAIL_REGULAR_EXPRESSION ".*[Tt]est: FAILED" 
  )
endforeach()

# source_group( dsxx\\test FILES ${test_lib_sources} ${test_sources} )


