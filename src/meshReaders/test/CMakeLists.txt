# ---------------------------------------------------------------------------- #
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
# ---------------------------------------------------------------------------- #
# $Id$
# ---------------------------------------------------------------------------- #

project( meshReaders_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_lib_sources 
   ${PROJECT_SOURCE_DIR}/meshReaders_test.cc 
   ${PROJECT_SOURCE_DIR}/meshReaders_test.hh
)
file( GLOB test_sources *.cc )
file( GLOB mesh_input_files *.mesh.in )
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

add_library( Lib_meshReaders_test ${DRACO_LIBRARY_TYPE} ${test_lib_sources} )
if( ${DRACO_LIBRARY_TYPE} MATCHES "SHARED" )
   set_target_properties( 
      Lib_meshReaders_test
        PROPERTIES 
        # Provide compile define macro to enable declspec(dllexport) linkage.
        COMPILE_DEFINITIONS BUILDING_DLL 
        # Use custom library naming
        OUTPUT ${CMAKE_SHARED_LIBRARY_PREFIX}meshReaders_test${CMAKE_SHARED_LIBRARY_SUFFIX}
   )
else()
   set_target_properties( 
      Lib_meshReaders_test 
      PROPERTIES 
        # Use custom library naming
	OUTPUT ${CMAKE_SHARED_LIBRARY_PREFIX}meshReaders_test${CMAKE_SHARED_LIBRARY_SUFFIX}
   )
endif()

# For Win32 with shared libraries, the package dll must be located in
# the test directory.

get_target_property( Lib_meshReaders_loc Lib_meshReaders LOCATION )
if( WIN32 )
   add_custom_command( TARGET Lib_meshReaders_test
      POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${Lib_meshReaders_loc} 
      ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
   )
endif()
add_dependencies( Lib_meshReaders_test Lib_meshReaders )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #
foreach( file ${test_sources} )
  get_filename_component( testname ${file} NAME_WE )
  add_executable( Ut_meshReaders_${testname}_exe ${file})
  set_target_properties( 
     Ut_meshReaders_${testname}_exe 
     PROPERTIES 
       OUTPUT_NAME ${testname} 
       VS_KEYWORD ${testname}
       PROJECT_LABEL Ut_meshReaders
  )
  target_link_libraries( Ut_meshReaders_${testname}_exe Lib_meshReaders Lib_meshReaders_test Lib_dsxx )
  add_test( meshReaders_${testname} ${testname} )
  set_tests_properties( meshReaders_${testname} PROPERTIES	
     PASS_REGULAR_EXPRESSION ".*[Tt]est: PASSED" 
     FAIL_REGULAR_EXPRESSION ".*[Tt]est: FAILED" 
  )
  if( "${test_loc}x" STREQUAL "x" )
    get_target_property( test_loc Ut_meshReaders_${testname}_exe LOCATION )
    get_filename_component( test_loc ${test_loc} PATH )
  endif()
endforeach()

# ---------------------------------------------------------------------------- #
# Copy mesh files to test directory (binary_dir)
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
add_custom_target( Ut_meshReaders_install_inputs ALL
    DEPENDS ${required_mesh_files}
)

