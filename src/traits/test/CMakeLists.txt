#=============================================================================
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
#=============================================================================
# $Id$
#=============================================================================

project( traits_test CXX )

# Find local headers
include_directories( 
   ${PROJECT_SOURCE_DIR}      # headers for tests
   ${PROJECT_SOURCE_DIR}/..   # headers for package
   ${MPI_INCLUDE_PATH} )

#set( test_lib_sources 
#	${PROJECT_SOURCE_DIR}/Traits_Test.hh)
file( GLOB test_sources *.cc )
#list( REMOVE_ITEM test_sources ${test_lib_sources} )

# Library Lib_traits_test
# add_library( Lib_traits_test ${DRACO_LIBRARY_TYPE} ${test_lib_sources} )
# if( ${DRACO_LIBRARY_TYPE} MATCHES "SHARED" )
	# set_target_properties( Lib_traits_test
		# PROPERTIES 
			# # Provide compile define macro to enable declspec(dllexport) linkage.
			# COMPILE_DEFINITIONS BUILDING_DLL 
			# # Use custom library naming
			# OUTPUT ${CMAKE_SHARED_LIBRARY_PREFIX}traits_test${CMAKE_SHARED_LIBRARY_SUFFIX}
	# )
# else()
	# set_target_properties( Lib_traits_test 
		# PROPERTIES 
			# # Use custom library naming
			# OUTPUT ${CMAKE_SHARED_LIBRARY_PREFIX}traits_test${CMAKE_SHARED_LIBRARY_SUFFIX}
	# )
# endif()
# For Win32 with shared libraries, the package dll must be located in the test directory.
# get_target_property( Lib_traits_loc Lib_traits LOCATION )
# if( WIN32 )
	# add_custom_command( TARGET Lib_traits_test
		# POST_BUILD
		# COMMAND ${CMAKE_COMMAND} -E copy_if_different ${Lib_traits_loc} 
		        # ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
		# )
# endif()
# add_dependencies( Lib_traits_test Lib_traits Lib_traits Lib_mesh_element Lib_c4 Lib_dsxx)

foreach( file ${test_sources} )
  get_filename_component( testname ${file} NAME_WE )
  add_executable( Ut_traits_${testname}_exe ${file})
  set_target_properties( Ut_traits_${testname}_exe 
   PROPERTIES 
      OUTPUT_NAME ${testname} 
      VS_KEYWORD ${testname}
      PROJECT_LABEL Ut_traits
      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  )
  target_link_libraries( Ut_traits_${testname}_exe 
    Lib_traits 
    #Lib_traits_test 
    Lib_dsxx )
endforeach()

foreach( file ${test_sources} )
  get_filename_component( testname ${file} NAME_WE )
  if( CMAKE_GENERATOR MATCHES "Visual Studio")
     set( test_loc  "${PROJECT_BINARY_DIR}/$(INTDIR)/${testname}" )
  else()
     get_target_property( test_loc Ut_traits_${testname}_exe LOCATION )
  endif()
  # Always run as scalar
  add_test( traits_${testname} ${testname} )
  set_tests_properties( traits_${testname} PROPERTIES	
     PASS_REGULAR_EXPRESSION ".*[Tt]est: PASSED"
     FAIL_REGULAR_EXPRESSION ".*[Tt]est: FAILED" 
     )
endforeach()

