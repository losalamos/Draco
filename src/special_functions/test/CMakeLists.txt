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
   ${ode_BINARY_DIR}          # ode/config.h
   ${GSL_INCLUDE_DIRS}
   )
# set( test_lib_sources 
file( GLOB test_sources *.cc )

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
      Lib_ode
      Lib_dsxx 
      ${GSL_LIBRARIES} )
   add_dependencies( Ut_special_functions_${testname}_exe Lib_special_functions)
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
