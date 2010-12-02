# ---------------------------------------------------------------------------- #
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
# ---------------------------------------------------------------------------- #
# $Id$
# ---------------------------------------------------------------------------- #

project( c4_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_lib_sources 
   ${PROJECT_SOURCE_DIR}/c4_test.cc 
   ${PROJECT_SOURCE_DIR}/c4_test.hh
)
file( GLOB test_sources *.cc )
list( REMOVE_ITEM test_sources ${test_lib_sources} )

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories( 
   ${PROJECT_SOURCE_DIR}      # headers for tests
   ${PROJECT_SOURCE_DIR}/..   # headers for package
   ${PROJECT_BINARY_DIR}/..   # config.h
   ${MPI_INCLUDE_PATH} )

# ---------------------------------------------------------------------------- #
# Build library for test directory
# ---------------------------------------------------------------------------- #

add_component_library( Lib_c4_test c4_test "${test_lib_sources}" )
add_dependencies( Lib_c4_test Lib_c4 )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

if( C4_MPI )
   set( NUM_PE_LIST 1;2;4 )
endif()

foreach( file ${test_sources} )
  get_filename_component( testname ${file} NAME_WE )
  add_executable( Ut_c4_${testname}_exe ${file})
  set_target_properties( 
    Ut_c4_${testname}_exe 
    PROPERTIES 
      OUTPUT_NAME ${testname} 
      VS_KEYWORD  ${testname}
      PROJECT_LABEL Ut_c4
      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  )
  target_link_libraries( Ut_c4_${testname}_exe 
     Lib_c4 
     Lib_c4_test 
     Lib_dsxx 
     ${MPI_LIBRARIES})
  # if( WIN32 )
   # add_custom_command( TARGET Ut_c4_${testname}_exe 
      # POST_BUILD
      # COMMAND ${CMAKE_COMMAND} -E copy_if_different  
      # ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
   # )
#endif()
endforeach()

if( C4_MPI )
   foreach( file ${test_sources} )
      get_filename_component( testname ${file} NAME_WE )
      #get_target_property( test_loc Ut_c4_${testname}_exe LOCATION )
      if( CMAKE_GENERATOR MATCHES "Visual Studio")
         set( test_loc  "${PROJECT_BINARY_DIR}/$(INTDIR)/${testname}" )
      else()
         get_target_property( test_loc Ut_c4_${testname}_exe LOCATION )
      endif()
      foreach( numPE ${NUM_PE_LIST} )
         if( "${testname}" MATCHES "phw" )
            add_test( 
               NAME c4_${testname}_${numPE} 
               COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${numPE} 
                       $<TARGET_FILE:Ut_c4_${testname}_exe> "hello" 
               )
            set_tests_properties( c4_${testname}_${numPE} PROPERTIES	
               PASS_REGULAR_EXPRESSION "Hello"
               FAIL_REGULAR_EXPRESSION "cruel" 
               )
         elseif( "${testname}" MATCHES "tstApplicationUnitTest" )
            add_test( 
                NAME c4_${testname}_${numPE} 
                COMMAND $<TARGET_FILE:Ut_c4_${testname}_exe> "--np" "${numPE}" )
            set_tests_properties( c4_${testname}_${numPE} PROPERTIES	
               PASS_REGULAR_EXPRESSION ".*[Tt]est: PASSED"
               FAIL_REGULAR_EXPRESSION ".*[Tt]est: FAILED" 
            )
         else()
            add_test( 
                NAME c4_${testname}_${numPE}
                COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${numPE} 
                        $<TARGET_FILE:Ut_c4_${testname}_exe> )
            set_tests_properties( c4_${testname}_${numPE} PROPERTIES	
               PASS_REGULAR_EXPRESSION ".*[Tt]est: PASSED"
               FAIL_REGULAR_EXPRESSION ".*[Tt]est: FAILED" 
            )
         endif()
        # add_dependencies( c4_${testname}_exe c4_test_lib )
      endforeach()
   endforeach()
else()
   # Scalar
   foreach( file ${test_sources} )
     if( "${testname}" MATCHES "phw" )
        add_test( c4_${testname} ${testname} "hello" )
        set_tests_properties( c4_${testname} PROPERTIES	
         PASS_REGULAR_EXPRESSION "Hello"
         FAIL_REGULAR_EXPRESSION "cruel" 
         )
     elseif( "${testname}" MATCHES "tstApplicationUnitTest" )
         message("${testname} scalr")
         add_test( c4_${testname}_${numPE} ${testname} "scalar" )
         set_tests_properties( c4_${testname}_${numPE} PROPERTIES	
            PASS_REGULAR_EXPRESSION ".*[Tt]est: PASSED"
            FAIL_REGULAR_EXPRESSION ".*[Tt]est: FAILED" 
         )
     else()
        add_test( c4_${testname} ${testname} )
        set_tests_properties( c4_${testname} PROPERTIES	
         PASS_REGULAR_EXPRESSION ".*[Tt]est: PASSED"
         FAIL_REGULAR_EXPRESSION ".*[Tt]est: FAILED" 
         )
      endif()
     # add_dependencies( c4_${testname}_exe c4_test_lib )
   endforeach()
endif()


   
