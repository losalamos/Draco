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

add_library( Lib_rng_test ${DRACO_LIBRARY_TYPE} ${test_lib_sources} )
if( ${DRACO_LIBRARY_TYPE} MATCHES "SHARED" )
   set_target_properties( 
      Lib_rng_test
        PROPERTIES 
        # Provide compile define macro to enable declspec(dllexport) linkage.
        COMPILE_DEFINITIONS BUILDING_DLL 
        # Use custom library naming
        OUTPUT ${CMAKE_SHARED_LIBRARY_PREFIX}rng_test${CMAKE_SHARED_LIBRARY_SUFFIX}
   )
else()
   set_target_properties( 
      Lib_rng_test 
      PROPERTIES 
        # Use custom library naming
	OUTPUT ${CMAKE_SHARED_LIBRARY_PREFIX}rng_test${CMAKE_SHARED_LIBRARY_SUFFIX}
   )
endif()

# For Win32 with shared libraries, the package dll must be located in
# the test directory.

get_target_property( Lib_rng_loc Lib_rng LOCATION )
if( WIN32 )
   add_custom_command( TARGET Lib_rng_test
      POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${Lib_rng_loc} 
      ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
   )
endif()
add_dependencies( Lib_rng_test Lib_rng )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

foreach( file ${test_sources} )
  get_filename_component( testname ${file} NAME_WE )
  add_executable( Ut_rng_${testname}_exe ${file})
  set_target_properties( 
     Ut_rng_${testname}_exe 
     PROPERTIES 
       OUTPUT_NAME ${testname} 
       VS_KEYWORD ${testname}
       PROJECT_LABEL Ut_rng
  )
  target_link_libraries( Ut_rng_${testname}_exe 
     Lib_rng 
     Lib_rng_test
     Lib_dsxx 
     ${GSL_LIBRARIES} )
  add_test( rng_${testname} ${testname} )
  set_tests_properties( rng_${testname} PROPERTIES	
     PASS_REGULAR_EXPRESSION ".*[Tt]est: PASSED" 
     FAIL_REGULAR_EXPRESSION ".*[Tt]est: FAILED" 
  )
endforeach()

# source_group( rng\\test FILES ${test_lib_sources} ${test_sources} )


