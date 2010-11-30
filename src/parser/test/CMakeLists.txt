#=============================================================================
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
#=============================================================================
# $Id$
#=============================================================================

project( parser_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB test_lib_sources *_pt.cc )
list( APPEND test_lib_sources 
   ${PROJECT_SOURCE_DIR}/parser_test.cc 
   ${PROJECT_SOURCE_DIR}/parser_test.hh)
file( GLOB test_sources *.cc )
list( REMOVE_ITEM test_sources ${test_lib_sources} )
file( GLOB input_files *.inp )

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories( 
   ${PROJECT_SOURCE_DIR}      # headers for tests
   ${PROJECT_SOURCE_DIR}/..   # headers for package
   ${PROJECT_BINARY_DIR}/..   # config.h
   ${MPI_INCLUDE_PATH} )

# Library Lib_parser_test
add_library( Lib_parser_test ${DRACO_LIBRARY_TYPE} ${test_lib_sources} )
if( ${DRACO_LIBRARY_TYPE} MATCHES "SHARED" )
   set_target_properties( Lib_parser_test
      PROPERTIES 
      # Provide compile define macro to enable declspec(dllexport) linkage.
      COMPILE_DEFINITIONS BUILDING_DLL 
      # Use custom library naming
      OUTPUT ${CMAKE_SHARED_LIBRARY_PREFIX}parser_test${CMAKE_SHARED_LIBRARY_SUFFIX}
      )
else()
   set_target_properties( Lib_parser_test 
      PROPERTIES 
      # Use custom library naming
      OUTPUT ${CMAKE_SHARED_LIBRARY_PREFIX}parser_test${CMAKE_SHARED_LIBRARY_SUFFIX}
      )
endif()
# For Win32 with shared libraries, the package dll must be located in the test directory.
get_target_property( Lib_parser_loc Lib_parser LOCATION )
if( WIN32 )
   add_custom_command( TARGET Lib_parser_test
      POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${Lib_parser_loc} 
      ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
      )
endif()
add_dependencies( Lib_parser_test Lib_parser )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

foreach( file ${test_sources} )
  get_filename_component( testname ${file} NAME_WE )
  add_executable( Ut_parser_${testname}_exe ${file})
  set_target_properties( Ut_parser_${testname}_exe 
   PROPERTIES 
      OUTPUT_NAME ${testname} 
      VS_KEYWORD ${testname}
      PROJECT_LABEL Ut_parser
      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  )
  target_link_libraries( Ut_parser_${testname}_exe 
    Lib_parser 
    Lib_parser_test 
    Lib_c4 
    Lib_dsxx 
    ${MPI_LIBRARIES} )
endforeach()

# ---------------------------------------------------------------------------- #
# Register Unit tests
# ---------------------------------------------------------------------------- #

if( ${C4_MPI} STREQUAL "1" )
   set( NUM_PE_LIST 1;2 )
endif()

if( ${C4_MPI} STREQUAL "1" )
   foreach( file ${test_sources} )
      get_filename_component( testname ${file} NAME_WE )
      if( CMAKE_GENERATOR MATCHES "Visual Studio")
         set( test_loc  "${PROJECT_BINARY_DIR}/$(INTDIR)/${testname}" )
      else()
         get_target_property( test_loc Ut_parser_${testname}_exe LOCATION )
      endif()
          # Always run as scalar
      if( "${testname}" MATCHES "tstToken" OR
          "${testname}" MATCHES "tstToken_Equivalence" OR
          "${testname}" MATCHES "tstutilities" OR
          "${testname}" MATCHES "tstUnit" OR
          "${testname}" MATCHES "tstString_Token_Stream" OR
          "${testname}" MATCHES "tstParse_Table" )
        add_test( parser_${testname} ${testname} )
        set_tests_properties( parser_${testname} PROPERTIES	
         PASS_REGULAR_EXPRESSION ".*[Tt]est: PASSED"
         FAIL_REGULAR_EXPRESSION ".*[Tt]est: FAILED" 
         )
      elseif( "${testname}" STREQUAL "tstConsole_Token_Stream" )
         # do nothing
      else() # Parallel tests
        foreach( numPE ${NUM_PE_LIST} )
          add_test( 
            NAME parser_${testname}_${numPE}
            COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${numPE} $<TARGET_FILE:Ut_parser_${testname}_exe> )
          set_tests_properties( parser_${testname}_${numPE} PROPERTIES	
            PASS_REGULAR_EXPRESSION ".*[Tt]est: PASSED"
            FAIL_REGULAR_EXPRESSION ".*[Tt]est: FAILED" 
          )
        endforeach()
      endif()
   endforeach()
else()
   # Scalar
   foreach( file ${test_sources} )
     if( "${testname}" MATCHES "tstToken" )
        add_test( parser_${testname} ${testname} "hello" )
        set_tests_properties( parser_${testname} PROPERTIES	
         PASS_REGULAR_EXPRESSION "Hello"
         FAIL_REGULAR_EXPRESSION "cruel" 
         )
     elseif( "${testname}" MATCHES "tstApplicationUnitTest" )
         add_test( parser_${testname}_${numPE} ${testname} "scalar" )
         set_tests_properties( parser_${testname}_${numPE} PROPERTIES	
            PASS_REGULAR_EXPRESSION ".*[Tt]est: PASSED"
            FAIL_REGULAR_EXPRESSION ".*[Tt]est: FAILED" 
         )
     else()
        add_test( parser_${testname} ${testname} )
        set_tests_properties( parser_${testname} PROPERTIES	
         PASS_REGULAR_EXPRESSION ".*[Tt]est: PASSED"
         FAIL_REGULAR_EXPRESSION ".*[Tt]est: FAILED" 
         )
      endif()
     # add_dependencies( parser_${testname}_exe parser_test_lib )
   endforeach()
endif()

# ---------------------------------------------------------------------------- #
# Copy input files to test directory (binary_dir)
# ---------------------------------------------------------------------------- #
unset(required_input_files)
foreach( inputfile ${input_files} )
    get_filename_component( inputfilenameonly ${inputfile} NAME_WE )
	if ( ${CMAKE_GENERATOR} MATCHES "Makefiles")
		set( outfile ${PROJECT_BINARY_DIR}/${inputfilenameonly}.inp )
	else ()
		set( outfile ${PROJECT_BINARY_DIR}/${CMAKE_BUILD_TYPE}/${inputfilenameonly}.inp )
	endif()
    add_custom_command( 
        OUTPUT  ${outfile}
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${inputfile} ${outfile}
        COMMENT "Copying ${inputfile} to ${outfile}"
    )
    list( APPEND required_input_files "${outfile}" )
endforeach()
add_custom_target( Ut_parser_install_inputs ALL
    DEPENDS ${required_input_files}
)
