# ---------------------------------------------------------------------------- #
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
# ---------------------------------------------------------------------------- #
# $Id$
# ---------------------------------------------------------------------------- #

project( linear_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_lib_sources 
   ${PROJECT_SOURCE_DIR}/linear_test.cc 
   ${PROJECT_SOURCE_DIR}/linear_test.hh
)
file( GLOB test_sources *.cc )
list( REMOVE_ITEM test_sources ${test_lib_sources} )
file( GLOB mesh_input_files *.in )
file( GLOB python_scripts python/*.py )
file( GLOB bench_files bench/*.0001 )

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

add_library( Lib_linear_test ${DRACO_LIBRARY_TYPE} ${test_lib_sources} )
if( ${DRACO_LIBRARY_TYPE} MATCHES "SHARED" )
   set_target_properties( 
      Lib_linear_test
        PROPERTIES 
        # Provide compile define macro to enable declspec(dllexport) linkage.
        COMPILE_DEFINITIONS BUILDING_DLL 
        # Use custom library naming
        OUTPUT ${CMAKE_SHARED_LIBRARY_PREFIX}linear_test${CMAKE_SHARED_LIBRARY_SUFFIX}
   )
else()
   set_target_properties( 
      Lib_linear_test 
      PROPERTIES 
        # Use custom library naming
	OUTPUT ${CMAKE_SHARED_LIBRARY_PREFIX}linear_test${CMAKE_SHARED_LIBRARY_SUFFIX}
   )
endif()

# For Win32 with shared libraries, the package dll must be located in
# the test directory.

get_target_property( Lib_linear_loc Lib_linear LOCATION )
if( WIN32 )
   add_custom_command( TARGET Lib_linear_test
      POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${Lib_linear_loc} 
      ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
   )
endif()
add_dependencies( Lib_linear_test Lib_linear )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

foreach( file ${test_sources} )
  get_filename_component( testname ${file} NAME_WE )
  add_executable( Ut_linear_${testname}_exe ${file})
  set_target_properties( 
     Ut_linear_${testname}_exe 
     PROPERTIES 
       OUTPUT_NAME ${testname} 
       VS_KEYWORD ${testname}
       PROJECT_LABEL Ut_linear
  )
  target_link_libraries( Ut_linear_${testname}_exe 
     Lib_linear 
     Lib_linear_test 
     Lib_traits 
     Lib_c4 
     Lib_dsxx 
     ${MPI_LIBRARY} )
  add_test( linear_${testname} ${testname} )
  set_tests_properties( linear_${testname} PROPERTIES	
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
    set( outmeshfile ${PROJECT_BINARY_DIR}/${meshfilenameonly} )
#    message("${meshfilenameonly} --> ${outmeshfile}")
    add_custom_command( 
        OUTPUT  ${outmeshfile}
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${meshfile} ${outmeshfile}
        COMMENT "Copying ${meshfile} to ${outmeshfile}"
    )
    list( APPEND required_mesh_files "${outmeshfile}" )
endforeach()
foreach( pythoninfile ${python_scripts} )
    get_filename_component( pythonfilenameonly ${pythoninfile} NAME_WE )
    set( pythonoutfile ${PROJECT_BINARY_DIR}/${pythonfilenameonly}.py )
    add_custom_command( 
        OUTPUT  ${pythonoutfile}
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${pythoninfile} ${pythonoutfile}
        COMMENT "Copying ${pythoninfile} to ${pythonoutfile}"
    )
    list( APPEND required_mesh_files "${pythonoutfile}" )
endforeach()
foreach( benchinfile ${bench_files} )
    get_filename_component( benchfilenameonly ${benchinfile} NAME )
    set( benchoutfile ${PROJECT_BINARY_DIR}/${benchfilenameonly} )
    add_custom_command( 
        OUTPUT  ${benchoutfile}
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${benchinfile} ${benchoutfile}
        COMMENT "Copying ${benchinfile} to ${benchoutfile}"
    )
    list( APPEND required_mesh_files "${benchoutfile}" )
endforeach()

add_custom_target( Ut_linear_install_inputs ALL
    DEPENDS ${required_mesh_files}
)



