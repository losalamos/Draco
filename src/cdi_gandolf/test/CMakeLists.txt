# ---------------------------------------------------------------------------- #
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
# ---------------------------------------------------------------------------- #
# $Id$
# ---------------------------------------------------------------------------- #

project( cdi_gandolf_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_lib_sources 
   ${PROJECT_SOURCE_DIR}/cdi_gandolf_test.cc 
   ${PROJECT_SOURCE_DIR}/cdi_gandolf_test.hh
)
file( GLOB test_sources *.cc )
list( REMOVE_ITEM test_sources ${test_lib_sources} )
file( GLOB ipcress_data *.ipcress )

##---------------------------------------------------------------------------##
## Copy data files
##---------------------------------------------------------------------------##

foreach( file ${ipcress_data} )
   get_filename_component( shortname ${file} NAME )
   configure_file( ${file} 
      ${PROJECT_BINARY_DIR}/${shortname} COPYONLY )
endforeach()

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

add_library( Lib_cdi_gandolf_test ${DRACO_LIBRARY_TYPE} ${test_lib_sources} )
if( ${DRACO_LIBRARY_TYPE} MATCHES "SHARED" )
   set_target_properties( 
      Lib_cdi_gandolf_test
        PROPERTIES 
        # Provide compile define macro to enable declspec(dllexport) linkage.
        COMPILE_DEFINITIONS BUILDING_DLL 
        # Use custom library naming
        OUTPUT ${CMAKE_SHARED_LIBRARY_PREFIX}cdi_gandolf_test${CMAKE_SHARED_LIBRARY_SUFFIX}
   )
else()
   set_target_properties( 
      Lib_cdi_gandolf_test 
      PROPERTIES 
        # Use custom library naming
	OUTPUT ${CMAKE_SHARED_LIBRARY_PREFIX}cdi_gandolf_test${CMAKE_SHARED_LIBRARY_SUFFIX}
   )
endif()

# For Win32 with shared libraries, the package dll must be located in
# the test directory.

get_target_property( Lib_cdi_gandolf_loc Lib_cdi_gandolf LOCATION )
if( WIN32 )
   add_custom_command( TARGET Lib_cdi_gandolf_test
      POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${Lib_cdi_gandolf_loc} 
      ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
   )
endif()
add_dependencies( Lib_cdi_gandolf_test Lib_cdi_gandolf )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

# Only add tests that are valid for this build configuration
if( rtt_cdi_gandolf_stub )
   set( test_sources ${PROJECT_SOURCE_DIR}/ReadOdfGandolfFile.cc )
endif()

foreach( file ${test_sources} )
  get_filename_component( testname ${file} NAME_WE )
  add_executable( Ut_cdi_gandolf_${testname}_exe ${file})
  set_target_properties( 
     Ut_cdi_gandolf_${testname}_exe 
     PROPERTIES 
       OUTPUT_NAME ${testname} 
       VS_KEYWORD ${testname}
       PROJECT_LABEL Ut_cdi_gandolf
  )
  target_link_libraries( Ut_cdi_gandolf_${testname}_exe 
    Lib_cdi_gandolf 
    Lib_cdi_gandolf_test 
    Lib_cdi 
    Lib_dsxx 
    ${GANDOLF_LIBRARIES}
  )
add_test( cdi_gandolf_${testname} ${testname} )
set_tests_properties( cdi_gandolf_${testname} PROPERTIES	
   PASS_REGULAR_EXPRESSION ".*[Tt]est: PASSED" 
   FAIL_REGULAR_EXPRESSION ".*[Tt]est: FAILED" 
   )
endforeach()

# source_group( cdi_gandolf\\test FILES ${test_lib_sources} ${test_sources} )


