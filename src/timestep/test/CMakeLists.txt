# ---------------------------------------------------------------------------- #
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
# ---------------------------------------------------------------------------- #
# $Id$
# ---------------------------------------------------------------------------- #

project( timestep_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_lib_sources 
   ${PROJECT_SOURCE_DIR}/timestep_test.cc 
   ${PROJECT_SOURCE_DIR}/timestep_test.hh
   ${PROJECT_SOURCE_DIR}/dummy_package.cc 
   ${PROJECT_SOURCE_DIR}/dummy_package.hh
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

add_library( Lib_timestep_test ${DRACO_LIBRARY_TYPE} ${test_lib_sources} )
if( ${DRACO_LIBRARY_TYPE} MATCHES "SHARED" )
   set_target_properties( 
      Lib_timestep_test
        PROPERTIES 
        # Provide compile define macro to enable declspec(dllexport) linkage.
        COMPILE_DEFINITIONS BUILDING_DLL 
        # Use custom library naming
        OUTPUT ${CMAKE_SHARED_LIBRARY_PREFIX}timestep_test${CMAKE_SHARED_LIBRARY_SUFFIX}
   )
else()
   set_target_properties( 
      Lib_timestep_test 
      PROPERTIES 
        # Use custom library naming
	OUTPUT ${CMAKE_SHARED_LIBRARY_PREFIX}timestep_test${CMAKE_SHARED_LIBRARY_SUFFIX}
   )
endif()

# For Win32 with shared libraries, the package dll must be located in
# the test directory.

get_target_property( Lib_timestep_loc Lib_timestep LOCATION )
if( WIN32 )
   add_custom_command( TARGET Lib_timestep_test
      POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${Lib_timestep_loc} 
      ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
   )
endif()
add_dependencies( Lib_timestep_test Lib_timestep )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

foreach( file ${test_sources} )
  get_filename_component( testname ${file} NAME_WE )
  add_executable( Ut_timestep_${testname}_exe ${file})
  set_target_properties( 
     Ut_timestep_${testname}_exe 
     PROPERTIES 
       OUTPUT_NAME ${testname} 
       VS_KEYWORD  ${testname}
       PROJECT_LABEL Ut_timestep
  )
  target_link_libraries( Ut_timestep_${testname}_exe 
  Lib_timestep_test
   Lib_timestep
   Lib_c4
   Lib_dsxx
   ${MPI_LIBRARIES} )
  # add_test( timestep_${testname} ${testname} )
  # set_tests_properties( timestep_${testname} PROPERTIES	
     # PASS_REGULAR_EXPRESSION ".*[Tt]est: PASSED" 
     # FAIL_REGULAR_EXPRESSION ".*[Tt]est: FAILED" 
  # )
endforeach()

# ---------------------------------------------------------------------------- #
# Register Unit tests
# ---------------------------------------------------------------------------- #

if( ${C4_MPI} STREQUAL "1" )
   set( NUM_PE_LIST 1;2;4 )
endif()

if( C4_MPI )
   foreach( file ${test_sources} )
      get_filename_component( testname ${file} NAME_WE )
      if( CMAKE_GENERATOR MATCHES "Visual Studio")
         set( test_loc  "${PROJECT_BINARY_DIR}/$(INTDIR)/${testname}" )
      else()
         get_target_property( test_loc Ut_parser_${testname}_exe LOCATION )
      endif()
      foreach( numPE ${NUM_PE_LIST} )
        add_test( 
          NAME    timestep_${testname}_${numPE}
          COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${numPE} 
                  $<TARGET_FILE:Ut_timestep_${testname}_exe> )
        set_tests_properties( timestep_${testname}_${numPE} PROPERTIES	
          PASS_REGULAR_EXPRESSION ".*[Tt]est: PASSED"
          FAIL_REGULAR_EXPRESSION ".*[Tt]est: FAILED" 
        )
      endforeach()
   endforeach()
else()
   # Scalar
   foreach( file ${test_sources} )
      add_test( parser_${testname} ${testname} )
      set_tests_properties( parser_${testname} PROPERTIES	
         PASS_REGULAR_EXPRESSION ".*[Tt]est: PASSED"
         FAIL_REGULAR_EXPRESSION ".*[Tt]est: FAILED" 
       )
   endforeach()
endif()
