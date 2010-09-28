#=============================================================================
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
#=============================================================================
# $Id$
#=============================================================================

project( special_functions_test CXX )

# Find local headers
include_directories( 
   ${PROJECT_SOURCE_DIR}      # headers for tests
   ${PROJECT_SOURCE_DIR}/..   # headers for package
   ${GSL_INCLUDE_DIRS}
)

# set( test_lib_sources 
	# ${PROJECT_SOURCE_DIR}/special_functions_test.cc 
	# ${PROJECT_SOURCE_DIR}/special_functions_test.hh)
file( GLOB test_sources *.cc )
# list( REMOVE_ITEM test_sources ${test_lib_sources} )
# file( GLOB input_files *.inp )

# Library Lib_special_functions_test
# add_library( Lib_special_functions_test ${DRACO_LIBRARY_TYPE} ${test_lib_sources} )
# if( ${DRACO_LIBRARY_TYPE} MATCHES "SHARED" )
	# set_target_properties( Lib_special_functions_test
		# PROPERTIES 
			# Provide compile define macro to enable declspec(dllexport) linkage.
			# COMPILE_DEFINITIONS BUILDING_DLL 
			# Use custom library naming
			# OUTPUT ${CMAKE_SHARED_LIBRARY_PREFIX}special_functions_test${CMAKE_SHARED_LIBRARY_SUFFIX}
	# )
# else()
	# set_target_properties( Lib_special_functions_test 
		# PROPERTIES 
			# Use custom library naming
			# OUTPUT ${CMAKE_SHARED_LIBRARY_PREFIX}special_functions_test${CMAKE_SHARED_LIBRARY_SUFFIX}
	# )
# endif()
# For Win32 with shared libraries, the package dll must be located in the test directory.
# get_target_property( Lib_special_functions_loc Lib_special_functions LOCATION )
# if( WIN32 )
	# add_custom_command( TARGET Lib_special_functions_test
		# POST_BUILD
		# COMMAND ${CMAKE_COMMAND} -E copy_if_different ${Lib_special_functions_loc} 
		        # ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
		# )
# endif()
# add_dependencies( Lib_special_functions_test Lib_special_functions )

foreach( file ${test_sources} )
  get_filename_component( testname ${file} NAME_WE )
  add_executable( Ut_special_functions_${testname}_exe ${file})
  set_target_properties( Ut_special_functions_${testname}_exe 
   PROPERTIES 
      OUTPUT_NAME ${testname} 
      VS_KEYWORD ${testname}
      PROJECT_LABEL Ut_special_functions
      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  )
  target_link_libraries( Ut_special_functions_${testname}_exe 
    Lib_special_functions 
    #Lib_special_functions_test 
    Lib_units
    Lib_dsxx 
    ${GSL_LIBRARIES} )
  add_test( special_functions_${testname} ${testname} )
  set_tests_properties( special_functions_${testname} PROPERTIES	
     PASS_REGULAR_EXPRESSION ".*[Tt]est: PASSED"
     FAIL_REGULAR_EXPRESSION ".*[Tt]est: FAILED" 
     )
endforeach()


# ---------------------------------------------------------------------------- #
# Copy vendor dll files to test directory (binary_dir)
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
add_custom_target( ut_special_functions_install_dll ALL
    DEPENDS ${required_input_files}
)
