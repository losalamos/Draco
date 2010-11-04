#=============================================================================
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
#=============================================================================
# $Id$
#=============================================================================

project( RTT_Format_Reader_test CXX )

# Find local headers
include_directories( 
   ${PROJECT_SOURCE_DIR}      # headers for tests
   ${PROJECT_SOURCE_DIR}/..   # headers for package
   ${PROJECT_BINARY_DIR}/..   # config.h
   ${MPI_INCLUDE_PATH} )

set( test_lib_sources 
	${PROJECT_SOURCE_DIR}/RTT_Format_Reader_test.cc 
	${PROJECT_SOURCE_DIR}/RTT_Format_Reader_test.hh)
file( GLOB test_sources *.cc )
list( REMOVE_ITEM test_sources ${test_lib_sources} )
file( GLOB mesh_input_files *.mesh.in )

# Library Lib_RTT_Format_Reader_test
add_library( Lib_RTT_Format_Reader_test ${DRACO_LIBRARY_TYPE} ${test_lib_sources} )
if( ${DRACO_LIBRARY_TYPE} MATCHES "SHARED" )
   set_target_properties( Lib_RTT_Format_Reader_test
      PROPERTIES 
      # Provide compile define macro to enable declspec(dllexport) linkage.
      COMPILE_DEFINITIONS BUILDING_DLL 
      # Use custom library naming
      OUTPUT ${CMAKE_SHARED_LIBRARY_PREFIX}RTT_Format_Reader_test${CMAKE_SHARED_LIBRARY_SUFFIX}
      )
else()
   set_target_properties( Lib_RTT_Format_Reader_test 
      PROPERTIES 
      # Use custom library naming
      OUTPUT ${CMAKE_SHARED_LIBRARY_PREFIX}RTT_Format_Reader_test${CMAKE_SHARED_LIBRARY_SUFFIX}
      )
endif()
# For Win32 with shared libraries, the package dll must be located in the test directory.
get_target_property( Lib_RTT_Format_Reader_loc Lib_RTT_Format_Reader LOCATION )
if( WIN32 )
   add_custom_command( TARGET Lib_RTT_Format_Reader_test
      POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${Lib_RTT_Format_Reader_loc} 
      ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
      )
endif()
add_dependencies( Lib_RTT_Format_Reader_test
   Lib_RTT_Format_Reader 
   Lib_RTT_Format_Reader
   Lib_mesh_element
   Lib_c4
   Lib_dsxx)

foreach( file ${test_sources} )
   get_filename_component( testname ${file} NAME_WE )
   add_executable( Ut_RTT_Format_Reader_${testname}_exe ${file})
  set_target_properties( Ut_RTT_Format_Reader_${testname}_exe 
   PROPERTIES 
      OUTPUT_NAME ${testname} 
      VS_KEYWORD ${testname}
      PROJECT_LABEL Ut_RTT_Format_Reader
      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  )
  target_link_libraries( Ut_RTT_Format_Reader_${testname}_exe 
    Lib_RTT_Format_Reader 
    Lib_RTT_Format_Reader_test 
    Lib_RTT_Format_Reader
    Lib_mesh_element
    Lib_c4 
    Lib_dsxx 
    ${MPI_LIBRARIES} )
endforeach()

foreach( file ${test_sources} )
   get_filename_component( testname ${file} NAME_WE )
   if( CMAKE_GENERATOR MATCHES "Visual Studio")
      set( test_loc  "${PROJECT_BINARY_DIR}/$(INTDIR)/${testname}" )
   else()
      get_target_property( test_loc Ut_RTT_Format_Reader_${testname}_exe LOCATION )
   endif()
   # Always run as scalar
   add_test( RTT_Format_Reader_${testname} ${testname} )
   set_tests_properties( RTT_Format_Reader_${testname} PROPERTIES	
      PASS_REGULAR_EXPRESSION ".*[Tt]est: PASSED"
      FAIL_REGULAR_EXPRESSION ".*[Tt]est: FAILED" 
      )
endforeach()


# ---------------------------------------------------------------------------- #
# Copy input files to test directory (binary_dir)
# ---------------------------------------------------------------------------- #
unset(required_mesh_files)
foreach( meshfile ${mesh_input_files} )
   get_filename_component( meshfilenameonly ${meshfile} NAME_WE )
   set( outmeshfile ${PROJECT_BINARY_DIR}/${meshfilenameonly}.mesh )
   add_custom_command( 
      OUTPUT  ${outmeshfile}
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${meshfile} ${outmeshfile}
      COMMENT "Copying ${meshfile} to ${outmeshfile}"
      )
   list( APPEND required_mesh_files "${outmeshfile}" )
endforeach()
add_custom_target( Ut_RTT_Format_Reader_install_inputs ALL
   DEPENDS ${required_mesh_files}
   )

