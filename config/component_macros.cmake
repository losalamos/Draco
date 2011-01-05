#-----------------------------*-cmake-*----------------------------------------#
# file   config/component_macros.cmake
# author 
# date   2010 Dec 1
# brief  Provide extra macros to simplify CMakeLists.txt for component
#        directories. 
# note   Copyright © 2010 LANS, LLC  
#------------------------------------------------------------------------------#
# $Id$ 
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------
# replacement for built in command 'add_library'
# 
# In addition to adding a library built from $sources, set
# Draco-specific properties for the library.  This macro reduces ~20
# lines of code down to 1-2.
# 
# Usage:
#
# add_component_library( <target_name> <output_library_name> "${list_of_sources}" )
#
# Note: you must use quotes around ${list_of_sources} to preserve the list.
#
# Example: see ds++/CMakeLists.txt
#
# Option: Consider using default_args (cmake.org/Wiki/CMakeMacroParseArguments)
#------------------------------------------------------------------------------
macro( add_component_library target_name outputname sources )

   add_library( ${target_name} ${DRACO_LIBRARY_TYPE} ${sources}  )
   if( "${DRACO_LIBRARY_TYPE}" MATCHES "SHARED" )
      set_target_properties( ${target_name} 
         PROPERTIES 
         # Provide compile define macro to enable declspec(dllexport) linkage.
         COMPILE_DEFINITIONS BUILDING_DLL 
         # Use custom library naming
         OUTPUT_NAME rtt_${outputname}
         )
   else()
      set_target_properties( ${target_name}
         PROPERTIES 
         # Use custom library naming
         OUTPUT_NAME rtt_${outputname}
         )
   endif()

   if( ${target_name} MATCHES "_test" )

      # This is a test library.  Find the component name
      string( REPLACE "_test" "" comp_target ${target_name} )

      # For Win32 with shared libraries, the package dll must be
      # located in the test directory.

      get_target_property( ${comp_target}_loc ${comp_target} LOCATION )
      if( WIN32 )
         add_custom_command( TARGET ${target_name}
            POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different ${${comp_target}_loc} 
                    ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
            )
      endif()

   endif()

   # OUTPUT_NAME ${CMAKE_STATIC_LIBRARY_PREFIX}rtt_ds++${CMAKE_STATIC_LIBRARY_SUFFIX}
   
endmacro()

# #----------------------------------------------------------------------#
# # add_ptest_contains
# #----------------------------------------------------------------------#
# macro( ptest_contains var value )
#    set( ${var} )
#    foreach( value2 ${ARGN} )
#       if( ${value} STREQUAL ${value2} )
#          set( ${var} TRUE )
#       endif()
#    endforeach()
# endmacro()

#----------------------------------------------------------------------#
# parse_arguments
# See cmake.org/Wiki/CMakeMacroParseArguments
#----------------------------------------------------------------------#

MACRO(PARSE_ARGUMENTS prefix arg_names option_names)
  SET(DEFAULT_ARGS)
  FOREACH(arg_name ${arg_names})    
    SET(${prefix}_${arg_name})
  ENDFOREACH(arg_name)
  FOREACH(option ${option_names})
    SET(${prefix}_${option} FALSE)
  ENDFOREACH(option)

  SET(current_arg_name DEFAULT_ARGS)
  SET(current_arg_list)
  FOREACH(arg ${ARGN})            
    SET(larg_names ${arg_names})    
    LIST(FIND larg_names "${arg}" is_arg_name)                   
    IF (is_arg_name GREATER -1)
      SET(${prefix}_${current_arg_name} ${current_arg_list})
      SET(current_arg_name ${arg})
      SET(current_arg_list)
    ELSE (is_arg_name GREATER -1)
      SET(loption_names ${option_names})    
      LIST(FIND loption_names "${arg}" is_option)            
      IF (is_option GREATER -1)
             SET(${prefix}_${arg} TRUE)
      ELSE (is_option GREATER -1)
             SET(current_arg_list ${current_arg_list} ${arg})
      ENDIF (is_option GREATER -1)
    ENDIF (is_arg_name GREATER -1)
  ENDFOREACH(arg)
  SET(${prefix}_${current_arg_name} ${current_arg_list})
ENDMACRO(PARSE_ARGUMENTS)


#----------------------------------------------------------------------#
# add_scalar_tests
#
# Given a list of sources, create unit test executables, one exe for
# each source file.  Register the test to be run by ctest.  
#
# Usage:
#
# add_scalar_tests( "${test_sources}" "${library_dependencies}" )
#
#----------------------------------------------------------------------#
macro( add_scalar_tests test_sources )

   parse_arguments( 
      # prefix
      addscalartest
      # list names
      "SOURCES;DEPS;TEST_ARGS;PASS_REGEX;FAIL_REGEX"
      # option names
      "NONE"
      ${ARGV}
      )

#    message("
# addscalartest_SOURCES    = ${addscalartest_SOURCES}
# addscalartest_DEPS       = ${addscalartest_DEPS}
# addscalartest_TEST_ARGS  = ${addscalartest_TEST_ARGS}
# addscalartest_PASS_REGEX = ${addscalartest_PASS_REGEX}
# addscalartest_FAIL_REGEX = ${addscalartest_FAIL_REGEX}
# ")

   # Sanity Checks
   if( "${addscalartest_SOURCES}none" STREQUAL "none" )
      message( FATAL_ERROR "You must provide the keyword SOURCES and a list of sources when using the add_scalar_tests macro.  Please see draco/config/component_macros.cmake::add_scalar_tests() for more information." )
   endif()

   # Pass/Fail criteria
   if( "${addscalartest_PASS_REGEX}none" STREQUAL "none" )
      set( addscalartest_PASS_REGEX ".*[Tt]est: PASSED" )
   endif()
   if( "${addscalartest_FAIL_REGEX}none" STREQUAL "none" )
      set( addscalartest_FAIL_REGEX ".*[Tt]est: FAILED" )
   endif()
 
   # What is the component name (always use Lib_${compname} as a dependency).
   string( REPLACE "_test" "" compname ${PROJECT_NAME} )
   set( iarg "0" )

   # If the test directory does not provide its own library (e.g.:
   # libc4_test.a), then don't try to link against it!
   get_target_property( test_lib_loc Lib_${compname}_test LOCATION )
   #message( "test_lib_loc = ${test_lib_loc}" )
   if( NOT "${test_lib_loc}" MATCHES "NOTFOUND" )
      set( test_lib_target_name "Lib_${compname}_test" )
      #message( "test_lib_target_name = ${test_lib_target_name}" )
   endif()

   # Loop over each test source files:
   # 1. Compile the executable
   # 2. Link against dependencies (libraries)
   # 3. Register the unit test
   # 4. Register the pass/fail criteria.

   foreach( file ${addscalartest_SOURCES} )
      get_filename_component( testname ${file} NAME_WE )
      add_executable( Ut_${compname}_${testname}_exe ${file} )
      set_target_properties( 
         Ut_${compname}_${testname}_exe 
         PROPERTIES 
           OUTPUT_NAME ${testname} 
           VS_KEYWORD  ${testname}
           PROJECT_LABEL Ut_${compname}
         )
      target_link_libraries( Ut_${compname}_${testname}_exe 
         # Lib_${compname}_test 
         ${test_lib_target_name}
         Lib_${compname} 
         ${addscalartest_DEPS}
         )
      if( "${addscalartest_TEST_ARGS}none" STREQUAL "none" )
         add_test( ${compname}_${testname} ${testname} )
         set_tests_properties( ${compname}_${testname} 
           PROPERTIES	
             PASS_REGULAR_EXPRESSION "${addscalartest_PASS_REGEX}"
             FAIL_REGULAR_EXPRESSION "${addscalartest_FAIL_REGEX}"
             )
       else()
          foreach( cmdarg ${addscalartest_TEST_ARGS} ) 
             math( EXPR iarg "${iarg} + 1" )
             separate_arguments( tmp UNIX_COMMAND ${cmdarg} )
             add_test( ${compname}_${testname}_arg${iarg} 
                ${testname} ${tmp} )
             # message("${compname}_${testname}_arg${iarg} ${testname} ${cmdarg}")
             set_tests_properties( ${compname}_${testname}_arg${iarg}
                PROPERTIES	
                  PASS_REGULAR_EXPRESSION "${addscalartest_PASS_REGEX}"
                  FAIL_REGULAR_EXPRESSION "${addscalartest_FAIL_REGEX}"
                )
          endforeach()
       endif()
   endforeach()

endmacro()

#----------------------------------------------------------------------#
# add_parallel_tests
#
# Given a list of sources, create unit test executables, one exe for
# each source file.  Register the test to be run by ctest.  
#
# Usage:
#
# add_scalar_tests( 
#    SOURCES "${test_sources}" 
#    DEPS    "${library_dependencies}" 
#    PE_LIST "1;2;4" )
#
#----------------------------------------------------------------------#
macro( add_parallel_tests )

   parse_arguments( 
      # prefix
      addparalleltest
      # list names
      "SOURCES;PE_LIST;DEPS;TEST_ARGS;PASS_REGEX;FAIL_REGEX"
      # option names
      "NONE"
      ${ARGV}
      )

   # Sanity Check
   if( "${addparalleltest_SOURCES}none" STREQUAL "none" )
      message( FATAL_ERROR "You must provide the keyword SOURCES and a list of sources when using the add_parallel_tests macro.  Please see draco/config/component_macros.cmake::add_parallel_tests() for more information." )
   endif()

   # Pass/Fail criteria
   if( "${addparalleltest_PASS_REGEX}none" STREQUAL "none" )
      set( addparalleltest_PASS_REGEX ".*[Tt]est: PASSED" )
   endif()
   if( "${addparalleltest_FAIL_REGEX}none" STREQUAL "none" )
      set( addparalleltest_FAIL_REGEX ".*[Tt]est: FAILED" )
   endif()

   # What is the component name (always use Lib_${compname} as a dependency).
   string( REPLACE "_test" "" compname ${PROJECT_NAME} )

   # Loop over each test source files:
   # 1. Compile the executable
   # 2. Link against dependencies (libraries)

   foreach( file ${addparalleltest_SOURCES} )
      get_filename_component( testname ${file} NAME_WE )
      add_executable( Ut_${compname}_${testname}_exe ${file} )
      set_target_properties( 
         Ut_${compname}_${testname}_exe 
         PROPERTIES 
           OUTPUT_NAME ${testname} 
           VS_KEYWORD  ${testname}
           PROJECT_LABEL Ut_${compname}
         )
      target_link_libraries( 
         Ut_${compname}_${testname}_exe 
         Lib_${compname}_test 
         Lib_${compname} 
         ${addparalleltest_DEPS}
         )
      # if( WIN32 )
      #    add_custom_command( TARGET Ut_c4_${testname}_exe 
      #       POST_BUILD
      #       COMMAND ${CMAKE_COMMAND} -E copy_if_different  
      #               ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
      #       )
      # endif()
   endforeach()

   # 3. Register the unit test
   # 4. Register the pass/fail criteria.
   if( C4_MPI )
      foreach( file ${addparalleltest_SOURCES} )
         get_filename_component( testname ${file} NAME_WE )
         if( CMAKE_GENERATOR MATCHES "Visual Studio")
            set( test_loc "${PROJECT_BINARY_DIR}/$(INTDIR)/${testname}" )
         else()
            get_target_property( test_loc 
               Ut_${compname}_${testname}_exe 
               LOCATION )
         endif()

         # Loop over PE_LIST, register test for each numPE
         foreach( numPE ${addparalleltest_PE_LIST} )
            add_test( 
               NAME    ${compname}_${testname}_${numPE}
               COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${numPE}
                       $<TARGET_FILE:Ut_${compname}_${testname}_exe> 
                       ${addparalleltest_TEST_ARGS}
               )
            set_tests_properties( ${compname}_${testname}_${numPE}
               PROPERTIES	
                 PASS_REGULAR_EXPRESSION "${addparalleltest_PASS_REGEX}"
                 FAIL_REGULAR_EXPRESSION "${addparalleltest_FAIL_REGEX}"
                 PROCESSORS              "${numPE}"
               )
         endforeach()
      endforeach()
   else( C4_MPI )
      # SCALAR Mode:
      foreach( file ${addparalleltest_SOURCES} )
         get_filename_component( testname ${file} NAME_WE )
         add_test( ${compname}_${testname} ${testname} 
                   ${addparalleltest_TEST_ARGS} )
         set_tests_properties( ${compname}_${testname} 
            PROPERTIES	
              PASS_REGULAR_EXPRESSION "${addparalleltest_PASS_REGEX}"
              FAIL_REGULAR_EXPRESSION "${addparalleltest_FAIL_REGEX}"
            )
      endforeach()
   endif( C4_MPI )

endmacro()

#----------------------------------------------------------------------#
#
#----------------------------------------------------------------------#
macro( provide_aux_files )

   parse_arguments( 
      # prefix
      auxfiles
      # list names
      "FILES;SRC_EXT;DEST_EXT"
      # option names
      "NONE"
      ${ARGV}
      )

   unset(required_files)
   foreach( file ${auxfiles_FILES} )
      get_filename_component( srcfilenameonly ${file} NAME )
      if( "${auxfiles_SRC_EXT}none" STREQUAL "none" )
         if( "${auxfiles_DEST_EXT}none" STREQUAL "none" )
            # do nothing
         else()
            # add DEST_EXT
            set( srcfilenameonly
               "${srcfilenameonly}${auxfiles_DEST_EXT}" )
         endif()
      else()
         if( "${auxfiles_DEST_EXT}none" STREQUAL "none" )
            # strip SRC_EXT
            string( REPLACE ${auxfiles_SRC_EXT} ""
               srcfilenameonly ${srcfilenameonly} )
         else()
            # replace SRC_EXT with DEST_EXT
            string( REPLACE ${auxfiles_SRC_EXT} ${auxfiles_DEST_EXT}
               srcfilenameonly ${srcfilenameonly} )
         endif()
      endif()
      if ( ${CMAKE_GENERATOR} MATCHES "Makefiles")
         set( outfile ${PROJECT_BINARY_DIR}/${srcfilenameonly} )
      else()
         set( outfile ${PROJECT_BINARY_DIR}/${CMAKE_BUILD_TYPE}/${srcfilenameonly} )
      endif()
      add_custom_command( 
         OUTPUT  ${outfile}
         COMMAND ${CMAKE_COMMAND} -E copy_if_different ${file} ${outfile}
         COMMENT "Copying ${file} to ${outfile}"
         )
      list( APPEND required_files "${outfile}" )
   endforeach()
   string( REPLACE "_test" "" compname ${PROJECT_NAME} )

   # Extra logic if multiple calls from the same directory.
   if( "${Ut_${compname}_install_inputs_iarg}notset" STREQUAL "notset" )
      set( Ut_${compname}_install_inputs_iarg "0" CACHE INTERNAL
   "counter for each provide_aux_files command.  Used to create individual targets for copying support files.")
   else()
      math( EXPR Ut_${compname}_install_inputs_iarg
         "${Ut_${compname}_install_inputs_iarg} + 1" )
   endif()
   add_custom_target(
      Ut_${compname}_install_inputs_${Ut_${compname}_install_inputs_iarg} 
      ALL
      DEPENDS ${required_files}
      )
   
endmacro()

#----------------------------------------------------------------------#
# End
#----------------------------------------------------------------------#
