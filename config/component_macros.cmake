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

   # This is a test library.  Find the component name
   string( REPLACE "_test" "" comp_target ${target_name} )
   # extract project name, minus leading "Lib_"
   string( REPLACE "Lib_" "" folder_name ${comp_target} )

   add_library( ${target_name} ${DRACO_LIBRARY_TYPE} ${sources}  )
   if( "${DRACO_LIBRARY_TYPE}" MATCHES "SHARED" )
      set_target_properties( ${target_name} 
         PROPERTIES 
         # Provide compile define macro to enable declspec(dllexport) linkage.
         COMPILE_DEFINITIONS BUILDING_DLL 
         # Use custom library naming
         OUTPUT_NAME rtt_${outputname}
         FOLDER ${folder_name}
         )
   else()
      set_target_properties( ${target_name}
         PROPERTIES 
         # Use custom library naming
         OUTPUT_NAME rtt_${outputname}
         FOLDER ${folder_name}
         )
   endif()

   if( ${target_name} MATCHES "_test" )
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

# ------------------------------------------------------------
# Register_scalar_test()
#
# 1. Special treatment for Roadrunner/ppe code (must ssh and then run)
# 2. Register the test
# 3. Register the pass/fail criteria.
# ------------------------------------------------------------
macro( register_scalar_test targetname runcmd command cmd_args )
   
   if( "${CMAKE_CXX_COMPILER}" MATCHES "ppu-g[+][+]" AND "${SITE}" MATCHES "rr[a-d][0-9]+a" )
      # Special treatment for Roadrunner PPE build.  The tests
      # must be run on the PPC chip on the backend.  If we are
      # running from the x86 backend, then we can run the tests
      # by ssh'ing to the 'b' node and running the test.
      add_test( 
         NAME    ${targetname}
         COMMAND ${RUN_CMD} "(cd ${PROJECT_BINARY_DIR};./${command} ${cmd_args})" )
   else()
      # Cielito needs the ./ in front of the binary name.
      if( "${MPIEXEC}" MATCHES "aprun" )
         set( APT_TARGET_FILE_PREFIX "./" )
      endif()
      separate_arguments( cmdargs UNIX_COMMAND ${cmd_args} )
      add_test( 
         NAME    ${targetname}
         COMMAND ${RUN_CMD} ${APT_TARGET_FILE_PREFIX}${command} ${cmdargs})
   endif()

   # set pass fail criteria, processors required, etc.
   set_tests_properties( ${targetname}
      PROPERTIES	
      PASS_REGULAR_EXPRESSION "${addscalartest_PASS_REGEX}"
      FAIL_REGULAR_EXPRESSION "${addscalartest_FAIL_REGEX}"
      PROCESSORS              "1"
      WORKING_DIRECTORY       "${PROJECT_BINARY_DIR}"
      )
   if( NOT "${addscalartest_RESOURCE_LOCK}none" STREQUAL "none" )
      set_tests_properties( ${targetname}
         PROPERTIES RESOURCE_LOCK "${addscalartest_RESOURCE_LOCK}" )
   endif()
   if( NOT "${addscalartest_RUN_AFTER}none" STREQUAL "none" )
      set_tests_properties( ${targetname}
         PROPERTIES DEPENDS "${addscalartest_RUN_AFTER}" )
   endif()
endmacro()

# ------------------------------------------------------------
# Register_parallel_test()
#
# 1. Register the test
# 2. Register the pass/fail criteria.
# ------------------------------------------------------------
macro( register_parallel_test targetname numPE command cmd_args )

   add_test( 
      NAME    ${targetname}
      COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${numPE}
              ${MPIRUN_POSTFLAGS}
              ${command}
              ${cmdarg}
      )
   set_tests_properties( ${targetname}
      PROPERTIES	
        PASS_REGULAR_EXPRESSION "${addparalleltest_PASS_REGEX}"
        FAIL_REGULAR_EXPRESSION "${addparalleltest_FAIL_REGEX}"
        PROCESSORS              "${numPE}"
        WORKING_DIRECTORY       "${PROJECT_BINARY_DIR}"
      )
   if( NOT "${addparalleltest_RESOURCE_LOCK}none" STREQUAL "none" )
      set_tests_properties( ${targetname}
         PROPERTIES RESOURCE_LOCK "${addparalleltest_RESOURCE_LOCK}" )
   endif()
   if( NOT "${addparalleltest_RUN_AFTER}none" STREQUAL "none" )
      set_tests_properties( ${targetname}
         PROPERTIES DEPENDS "${addparalleltest_RUN_AFTER}" )
   endif()

endmacro()

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

   # These become variables of the form ${addscalartests_SOURCES}, etc.
   parse_arguments( 
      # prefix
      addscalartest
      # list names
      "SOURCES;DEPS;TEST_ARGS;PASS_REGEX;FAIL_REGEX;RESOURCE_LOCK;RUN_AFTER"
      # option names
      "NONE"
      ${ARGV}
      )

   # Special Cases:
   # ------------------------------------------------------------
   # On some platforms (Cielo, Cielito, RedStorm), even scalar tests
   # must be run underneath MPIEXEC (yod, aprun):
   if( "${MPIEXEC}" MATCHES "aprun" )
      set( RUN_CMD ${MPIEXEC} -n 1 )
      set( APT_TARGET_FILE_PREFIX "./" )
   elseif(  "${MPIEXEC}" MATCHES "yod" )
      set( RUN_CMD ${MPIEXEC} -np 1 )
   else()
      unset( RUN_CMD )
   endif()

   # When on roadrunner backend (x86), we must ssh into the PPE to run
   # the test.  If we are on roadrunner frontend, we cannot run the
   # tests.
   if( "${CMAKE_CXX_COMPILER}" MATCHES "ppu-g[+][+]" )
      if( "${SITE}" MATCHES "rr[a-d][0-9]+a" )
         string( REGEX REPLACE "a[.]rr[.]lanl[.]gov" "b" ppe_node ${SITE} )
         set( RUN_CMD ssh ${ppe_node} )
      endif()
   endif()

   # Sanity Checks
   # ------------------------------------------------------------
   if( "${addscalartest_SOURCES}none" STREQUAL "none" )
      message( FATAL_ERROR "You must provide the keyword SOURCES and a list of sources when using the add_scalar_tests macro.  Please see draco/config/component_macros.cmake::add_scalar_tests() for more information." )
   endif()

   # Pass/Fail criteria
   if( "${addscalartest_PASS_REGEX}none" STREQUAL "none" )
      set( addscalartest_PASS_REGEX ".*[Tt]est: PASSED" )
   endif()
   if( "${addscalartest_FAIL_REGEX}none" STREQUAL "none" )
      set( addscalartest_FAIL_REGEX ".*[Tt]est: FAILED" )
      list( APPEND addscalartest_FAIL_REGEX ".*ERROR:.*" )
   endif()
 
   # Format resource lock command
   if( NOT "${addscalartest_RESOURCE_LOCK}none" STREQUAL "none" )
      set( addscalartest_RESOURCE_LOCK 
         "RESOURCE_LOCK ${addscalartest_RESOURCE_LOCK}")
   endif()

   # What is the component name (always use Lib_${compname} as a dependency).
   string( REPLACE "_test" "" compname ${PROJECT_NAME} )
   set( iarg "0" )

   # If the test directory does not provide its own library (e.g.:
   # libc4_test.a), then don't try to link against it!
   get_target_property( test_lib_loc Lib_${compname}_test LOCATION )
   if( NOT "${test_lib_loc}" MATCHES "NOTFOUND" )
      set( test_lib_target_name "Lib_${compname}_test" )
   endif()
   get_target_property( pkg_lib_loc Lib_${compname} LOCATION )
   if( NOT "${pkg_lib_loc}" MATCHES "NOTFOUND" )
      list( APPEND test_lib_target_name "Lib_${compname}" )
   endif()

   # Loop over each test source files:
   # 1. Compile the executable
   # 3. Register the unit test

   # Generate the executable
   # ------------------------------------------------------------
   foreach( file ${addscalartest_SOURCES} )

      if( "${file}" MATCHES "tstParallelUnitTest" )
         message("RUN_CMD = ${RUN_CMD}")
      endif()

      get_filename_component( testname ${file} NAME_WE )
      add_executable( Ut_${compname}_${testname}_exe ${file} )
      set_target_properties( 
         Ut_${compname}_${testname}_exe 
         PROPERTIES 
           OUTPUT_NAME ${testname} 
           VS_KEYWORD  ${testname}
           FOLDER ${compname}
         )
      target_link_libraries( 
         Ut_${compname}_${testname}_exe 
         ${test_lib_target_name}
         ${addscalartest_DEPS}
         )
   endforeach()

   # Register the unit test
   # ------------------------------------------------------------
   foreach( file ${addscalartest_SOURCES} )
      get_filename_component( testname ${file} NAME_WE )

      if( "${addscalartest_TEST_ARGS}none" STREQUAL "none" )
         register_scalar_test( ${compname}_${testname} 
            "${RUN_CMD}" ${testname} "" )
       else()
          foreach( cmdarg ${addscalartest_TEST_ARGS} ) 
             math( EXPR iarg "${iarg} + 1" )
             register_scalar_test( ${compname}_${testname}_arg${iarg} 
                "${RUN_CMD}" ${testname} "${cmdarg}" )
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
# add_parallel_tests( 
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
      "SOURCES;PE_LIST;DEPS;TEST_ARGS;PASS_REGEX;FAIL_REGEX;RESOURCE_LOCK;RUN_AFTER;MPIFLAGS"
      # option names
      "NONE"
      ${ARGV}
      )

   # When on roadrunner backend (x86), we must ssh into the PPE to run
   # the test.  If we are on roadrunner frontend, we cannot run the
   # tests.
   if( "${CMAKE_CXX_COMPILER}" MATCHES "ppu-g[+][+]" )
      if( "${SITE}" MATCHES "rr[a-d][0-9]+a" )
         string( REGEX REPLACE "a[.]rr[.]lanl[.]gov" "b" ppe_node ${SITE} )
         set( RUN_CMD ssh ${ppe_node} )
      endif()
   endif()

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
      list( APPEND addparalleltest_FAIL_REGEX ".*ERROR:.*" )
   endif()

   # Format resource lock command
   if( NOT "${addparalleltest_RESOURCE_LOCK}none" STREQUAL "none" )
      set( addparalleltest_RESOURCE_LOCK 
         "RESOURCE_LOCK ${addparalleltest_RESOURCE_LOCK}")
   endif()

   # What is the component name (always use Lib_${compname} as a dependency).
   string( REPLACE "_test" "" compname ${PROJECT_NAME} )

   # If the test directory does not provide its own library (e.g.:
   # libc4_test.a), then don't try to link against it!
   get_target_property( test_lib_loc Lib_${compname}_test LOCATION )
   if( NOT "${test_lib_loc}" MATCHES "NOTFOUND" )
      set( test_lib_target_name "Lib_${compname}_test" )
   endif()

   # Override MPI Flags upon user request
   if ( NOT DEFINED addparalleltest_MPIFLAGS )
       set( MPIRUN_POSTFLAGS ${MPIEXEC_POSTFLAGS} )
   else()
       set( MPIRUN_POSTFLAGS "${addparalleltest_MPIFLAGS}" )
   endif ()


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
           FOLDER ${compname}
         )
      target_link_libraries( 
         Ut_${compname}_${testname}_exe 
         ${test_lib_target_name} 
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
   if( ${DRACO_C4} MATCHES "MPI" )
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

         # 2011-02-22 KT: I noticed that a handful of tests were
         # failing because MPI was aborting.  After some investigation
         # this appears be be a result of how OpenMPI is installed on
         # the big iron.  That is, our installations do not allow
         # multiple mpirun executing on the same core. 
         #
         # In an attempt to fix this problem, I added "--bind-to-none"
         # to the add_test() command.  This does allow multiple mpirun
         # on the same core.  However, each mpirun job was using the
         # same 4 cores while the remaining 12 where idle.  
         #
         # Jon Dahl pointed me to the "--mca mpi_paffinity_alone 0"
         # option that appears to allow me to run on all 16 cores and
         # allows multiple mpirun to execute at the same time on each
         # core.  This appears to fix the randomly failing tests.
         #
         # http://www.open-mpi.org/faq/?category=tuning#using-paffinity-v1.4

         foreach( numPE ${addparalleltest_PE_LIST} )
            set( iarg 0 )
            if( "${addparalleltest_TEST_ARGS}none" STREQUAL "none" )
               register_parallel_test( 
                  ${compname}_${testname}_${numPE}
                  ${numPE}
                  $<TARGET_FILE:Ut_${compname}_${testname}_exe>
                  "" )
            else()
               foreach( cmdarg ${addparalleltest_TEST_ARGS} ) 
                  math( EXPR iarg "${iarg} + 1" )
                  register_parallel_test( 
                     ${compname}_${testname}_${numPE}_arg${iarg}
                     ${numPE} 
                     $<TARGET_FILE:Ut_${compname}_${testname}_exe>
                     ${cmdarg} )
               endforeach()
            endif()
         endforeach()
      endforeach()
   else( ${DRACO_C4} MATCHES "MPI" )
      # SCALAR Mode:
      foreach( file ${addparalleltest_SOURCES} )
         set( iarg "0" )
         get_filename_component( testname ${file} NAME_WE )

         set( addscalartest_PASS_REGEX "${addparalleltest_PASS_REGEX}" )
         set( addscalartest_FAIL_REGEX "${addparalleltest_FAIL_REGEX}" )
         set( addscalartest_RESOURCE_LOCK "${addparalleltest_RESOURCE_LOCK}" )
         set( addscalartest_RUN_AFTER "${addparalleltest_RUN_AFTER}" )

         if( "${addparalleltest_TEST_ARGS}none" STREQUAL "none" )
            
            register_scalar_test( ${compname}_${testname} 
               "${RUN_CMD}" ${testname} "" )
         else()

            foreach( cmdarg ${addparalleltest_TEST_ARGS} ) 
               math( EXPR iarg "${iarg} + 1" )
               register_scalar_test( ${compname}_${testname}_arg${iarg}
                  "${RUN_CMD}" ${testname} "${cmdarg}" )
            endforeach()

         endif()

      endforeach()
   endif( ${DRACO_C4} MATCHES "MPI" )

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

   # Strip any CVS directories from aux. files list
   foreach( file ${auxfiles_FILES} )
      if( ${file} MATCHES CVS$ )
         list( REMOVE_ITEM auxfiles_FILES ${file})
      endif()
   endforeach()
 

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
      if ( ${CMAKE_GENERATOR} MATCHES "Makefiles" OR 
           ${CMAKE_GENERATOR} MATCHES "Xcode"  )
         set( outfile ${PROJECT_BINARY_DIR}/${srcfilenameonly} )
      else()
         set( outfile ${PROJECT_BINARY_DIR}/${CMAKE_BUILD_TYPE}/${srcfilenameonly} )
      endif()
      add_custom_command( 
         OUTPUT  ${outfile}
         COMMAND ${CMAKE_COMMAND} -E copy_if_different ${file} ${outfile}
         DEPENDS ${file}
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
# CONDITIONALLY_ADD_SUBDIRECTORY - add a directory to the build while 
#      allowing exceptions:
# 
# E.g.: conditionally_add_subdirectory( 
#      COMPONENT "mc"
#      CXX_COMPILER_EXCEPTION "spu-g[+][+]" )
#----------------------------------------------------------------------#
macro( conditionally_add_subdirectory )

   parse_arguments( 
      # prefix
      caddsubdir
      # list names
      "COMPONENTS;CXX_COMPILER_EXCEPTION;CXX_COMPILER_MATCHES"
      # option names
      "NONE"
      ${ARGV}
      )

   # if the current compiler doesn't match the provided regex, then
   # add the directory to the build
   if( caddsubdir_CXX_COMPILER_EXCEPTION )
      if( NOT "${CMAKE_CXX_COMPILER}" MATCHES
            "${caddsubdir_CXX_COMPILER_EXCEPTION}" )
         foreach( comp ${caddsubdir_COMPONENTS} )
            message(STATUS "Configuring ${comp}")
            add_subdirectory( ${comp} )
         endforeach()
      endif()
   endif()

   # Only add the component if the current compiler matches the
   # requested regex
   if( caddsubdir_CXX_COMPILER_MATCHES )
      if( "${CMAKE_CXX_COMPILER}" MATCHES
            "${caddsubdir_CXX_COMPILER_MATCHES}" )
         foreach( comp ${caddsubdir_COMPONENTS} )
            message(STATUS "Configuring ${comp}")
            add_subdirectory( ${comp} )
         endforeach()
      endif()
   endif()


endmacro()

#----------------------------------------------------------------------#
# End
#----------------------------------------------------------------------#
