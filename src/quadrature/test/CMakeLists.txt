#=============================================================================
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
#=============================================================================
# $Id$
#=============================================================================

project( quadrature_test CXX )

# Find local headers
include_directories( 
   ${PROJECT_SOURCE_DIR}      # headers for tests
   ${PROJECT_SOURCE_DIR}/..   # headers for package
   ${PROJECT_BINARY_DIR}/..   # config.h
   ${MPI_INCLUDE_PATH} )

set( test_lib_sources 
	${PROJECT_SOURCE_DIR}/quadrature_test.cc 
	${PROJECT_SOURCE_DIR}/quadrature_test.hh)
file( GLOB test_sources *.cc )
list( REMOVE_ITEM test_sources ${test_lib_sources} )
file( GLOB input_files *.inp )

# Library Lib_quadrature_test
add_library( Lib_quadrature_test ${DRACO_LIBRARY_TYPE} ${test_lib_sources} )
if( ${DRACO_LIBRARY_TYPE} MATCHES "SHARED" )
   set_target_properties( Lib_quadrature_test
      PROPERTIES 
        # Provide compile define macro to enable declspec(dllexport) linkage.
        COMPILE_DEFINITIONS BUILDING_DLL 
        # Use custom library naming
        OUTPUT ${CMAKE_SHARED_LIBRARY_PREFIX}quadrature_test${CMAKE_SHARED_LIBRARY_SUFFIX}
      )
else()
   set_target_properties( Lib_quadrature_test 
      PROPERTIES 
        # Use custom library naming
        OUTPUT ${CMAKE_SHARED_LIBRARY_PREFIX}quadrature_test${CMAKE_SHARED_LIBRARY_SUFFIX}
      )
endif()
# For Win32 with shared libraries, the package dll must be located in the test directory.
get_target_property( Lib_quadrature_loc Lib_quadrature LOCATION )
if( WIN32 )
   add_custom_command( TARGET Lib_quadrature_test
      POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${Lib_quadrature_loc} 
      ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
      )
endif()
add_dependencies( Lib_quadrature_test Lib_quadrature )

foreach( file ${test_sources} )
  get_filename_component( testname ${file} NAME_WE )
  add_executable( Ut_quadrature_${testname}_exe ${file})
  set_target_properties( Ut_quadrature_${testname}_exe 
   PROPERTIES 
      OUTPUT_NAME ${testname} 
      VS_KEYWORD ${testname}
      PROJECT_LABEL Ut_quadrature
      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  )
  target_link_libraries( Ut_quadrature_${testname}_exe 
    Lib_quadrature 
    Lib_quadrature_test 
    Lib_parser
    Lib_c4 
    Lib_special_functions
    Lib_dsxx 
    ${GSL_LIBRARIES}
    ${MPI_LIBRARIES} )
endforeach()

foreach( file ${test_sources} )
  get_filename_component( testname ${file} NAME_WE )
  if( CMAKE_GENERATOR MATCHES "Visual Studio")
     set( test_loc  "${PROJECT_BINARY_DIR}/$(INTDIR)/${testname}" )
  else()
     get_target_property( test_loc Ut_quadrature_${testname}_exe LOCATION )
  endif()
  # Always run as scalar
  add_test( quadrature_${testname} ${testname} )
  set_tests_properties( quadrature_${testname} PROPERTIES	
     PASS_REGULAR_EXPRESSION ".*[Tt]est: PASSED"
     FAIL_REGULAR_EXPRESSION ".*[Tt]est: FAILED" 
     )
endforeach()


# ---------------------------------------------------------------------------- #
# Copy input files to test directory (binary_dir)
# ---------------------------------------------------------------------------- #
unset(required_input_files)
foreach( dllfile ${GSL_DLL_LIBRARIES} )
    get_filename_component( dllfilenameonly ${dllfile} NAME_WE )
    set( outfile ${PROJECT_BINARY_DIR}/${dllfilenameonly}.dll )
    add_custom_command( 
        OUTPUT  ${outfile}
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${dllfile} ${outfile}
        COMMENT "Copying ${dllfile} to ${outfile}"
    )
    list( APPEND required_input_files "${outfile}" )
endforeach()
add_custom_target( ut_quadrature_install_dll ALL
    DEPENDS ${required_input_files}
)

