#-----------------------------*-cmake-*----------------------------------------#
# file   config/component_macros.cmake
# author Kelly G. Thompson, kgt@lanl.gov
# date   2010 Dec 1
# brief  Provide extra macros to simplify CMakeLists.txt for component
#        directories. 
# note   Copyright (C) 2010-2014 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$ 
#------------------------------------------------------------------------------#

# requires parse_arguments()
include( parse_arguments )

#------------------------------------------------------------------------------
# replacement for built in command 'add_library'
# 
# Purpose 1: In addition to adding a library built from $sources, set
# Draco-specific properties for the library.  This macro reduces ~20
# lines of code down to 1-2.
#
# Purpose 2: Encapsulate library and vendor library dependencies per
# package.
#
# Purpose 3: Use information from 1 and 2 above to generate exported
# targets. 
# 
# Usage:
#
# add_component_library(
#   TARGET       "target name"
#   LIBRARY_NAME "output library name"
#   TARGET_DEPS  "dep1;dep2;..."
#   PREFIX       "ClubIMC"
#   SOURCES      "file1.cc;file2.cc;..."
#   LIBRARY_NAME_PREFIX "rtt_"
#   VENDOR_LIST  "MPI;GSL"
#   VENDOR_LIBS  "${MPI_CXX_LIBRARIES};${GSL_LIBRARIES}"
#   VENDOR_INCLUDE_DIRS "${MPI_CXX_INCLUDE_DIR};${GSL_INCLUDE_DIR}"
#   )
#
# Example:
#
# add_component_library(
#   TARGET       Lib_quadrature
#   LIBRARY_NAME quadrature
#   TARGET_DEPS  "Lib_parser;Lib_special_functions;Lib_mesh_element"
#   PREFIX       "Draco"
#   SOURCES      "${sources}"
#   )
#
# Note: you must use quotes around ${list_of_sources} to preserve the list.
#------------------------------------------------------------------------------
macro( add_component_library )
#target_name outputname sources 
# optional argument: libraryPrefix 

 # These become variables of the form ${acl_NAME}, etc.
   parse_arguments( 
      # prefix
      acl
      # list names
      "PREFIX;TARGET;LIBRARY_NAME;SOURCES;TARGET_DEPS;VENDOR_LIST;VENDOR_LIBS;VENDOR_INCLUDE_DIRS;LIBRARY_NAME_PREFIX;LINK_LANGUAGE"
      # option names
      "NOEXPORT"
      ${ARGV}
      )

   #
   # Defaults:
   # 
   # Optional 3rd argument is the library prefix.  The default is "rtt_".
   if( "${acl_LIBRARY_NAME_PREFIX}x" STREQUAL "x" )
      set( acl_LIBRARY_NAME_PREFIX "rtt_" )
   endif()
   if( "${acl_LINK_LANGUAGE}x" STREQUAL "x" )
      set( acl_LINK_LANGUAGE CXX )
   endif()

   #
   # Create the Library and set the Properties
   #

   # This is a test library.  Find the component name
   string( REPLACE "_test" "" comp_target ${acl_TARGET} )
   # extract project name, minus leading "Lib_"
   string( REPLACE "Lib_" "" folder_name ${acl_TARGET} )

   add_library( ${acl_TARGET} ${DRACO_LIBRARY_TYPE} ${acl_SOURCES} )
   if( "${DRACO_LIBRARY_TYPE}" MATCHES "SHARED" )
      set( compdefs COMPILE_DEFINITIONS BUILDING_DLL )
   endif()
   set_target_properties( ${acl_TARGET} 
      PROPERTIES 
      # Provide compile define macro to enable declspec(dllexport) linkage.
      ${compdefs}
      # COMPILE_DEFINITIONS BUILDING_DLL
      # Use custom library naming
      OUTPUT_NAME ${acl_LIBRARY_NAME_PREFIX}${acl_LIBRARY_NAME}
      FOLDER      ${folder_name}
      )

   #
   # Generate properties related to library dependencies
   #
   if( NOT "${acl_TARGET_DEPS}x" STREQUAL "x" )
      target_link_libraries( ${acl_TARGET} ${acl_TARGET_DEPS} )
   endif()
   if( NOT "${acl_VENDOR_LIBS}x" STREQUAL "x" )
      target_link_libraries( ${acl_TARGET} ${acl_VENDOR_LIBS} )
   endif()
   if( NOT "${acl_VENDOR_INCLUDE_DIRS}x" STREQUAL "x" )
      include_directories( ${acl_VENDOR_INCLUDE_DIRS} )
   endif()

   #
   # Register the library for exported library support
   #
   
   # Defaults
   if( "${acl_PREFIX}x" STREQUAL "x" )
      set( acl_PREFIX "Draco" )
   endif()

   # Find target file name and location
   get_target_property( impname ${acl_TARGET} OUTPUT_NAME )

   # the above command returns the location in the build tree.  We
   # need to convert this to the install location.
   # get_filename_component( imploc ${imploc} NAME )
   if( ${DRACO_SHARED_LIBS} )
      set( imploc "${CMAKE_INSTALL_PREFIX}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}${impname}${CMAKE_SHARED_LIBRARY_SUFFIX}" )
   else()
      set( imploc "${CMAKE_INSTALL_PREFIX}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}${impname}${CMAKE_STATIC_LIBRARY_SUFFIX}" )
   endif()

   set( ilil "")
   if( "${acl_TARGET_DEPS}x" STREQUAL "x" AND  "${acl_VENDOR_LIBS}x" STREQUAL "x")
      # do nothing
   elseif( "${acl_TARGET_DEPS}x" STREQUAL "x" )
      set( ilil "${acl_VENDOR_LIBS}" )
   elseif( "${acl_VENDOR_LIBS}x" STREQUAL "x")
      set( ilil "${acl_TARGET_DEPS}" )
   else()
      set( ilil "${acl_TARGET_DEPS};${acl_VENDOR_LIBS}" )
   endif()
   
   # For non-test libraries, save properties to the
   # project-config.cmake file     
   if( "${ilil}x" STREQUAL "x" )
     set( ${acl_PREFIX}_EXPORT_TARGET_PROPERTIES 
         "${${acl_PREFIX}_EXPORT_TARGET_PROPERTIES}
set_target_properties(${acl_TARGET} PROPERTIES
   IMPORTED_LINK_INTERFACE_LANGUAGES \"${acl_LINK_LANGUAGE}\"
   IMPORTED_LOCATION                 \"${imploc}\"
)
")
   else()
     set( ${acl_PREFIX}_EXPORT_TARGET_PROPERTIES 
       "${${acl_PREFIX}_EXPORT_TARGET_PROPERTIES}
set_target_properties(${acl_TARGET} PROPERTIES
   IMPORTED_LINK_INTERFACE_LANGUAGES \"${acl_LINK_LANGUAGE}\"
   IMPORTED_LINK_INTERFACE_LIBRARIES \"${ilil}\"
   IMPORTED_LOCATION                 \"${imploc}\"
)
")
   endif()
 
   # Set the install rpath for apple builds
   # if( ${DRACO_SHARED_LIBS} )
   #   set_target_properties( ${acl_TARGET} PROPERTIES INSTALL_RPATH "${imploc}" )
   # endif()

   # Only publish information to draco-config.cmake for non-test
   # libraries.  Also, omit any libraries that are marked as NOEXPORT
   if( NOT ${acl_NOEXPORT} AND
       NOT "${acl_TARGET}" MATCHES "test" )

      list( APPEND ${acl_PREFIX}_LIBRARIES ${acl_TARGET} )
      string( REPLACE "Lib_" "" compname ${acl_TARGET} )
      list( APPEND ${acl_PREFIX}_PACKAGE_LIST ${compname} )

      list( APPEND ${acl_PREFIX}_TPL_LIST ${acl_VENDOR_LIST} )
      list( APPEND ${acl_PREFIX}_TPL_INCLUDE_DIRS ${acl_VENDOR_INCLUDE_DIRS} )
      list( APPEND ${acl_PREFIX}_TPL_LIBRARIES ${acl_VENDOR_LIBS} )
      if( ${acl_PREFIX}_TPL_INCLUDE_DIRS )
         list( REMOVE_DUPLICATES ${acl_PREFIX}_TPL_INCLUDE_DIRS )
      endif()
      if( ${acl_PREFIX}_TPL_LIBRARIES )
         list( REMOVE_DUPLICATES ${acl_PREFIX}_TPL_LIBRARIES )
      endif()
      if( ${acl_PREFIX}_TPL_LIST )
         list( REMOVE_DUPLICATES ${acl_PREFIX}_TPL_LIST )   
      endif()
      
      set( ${acl_PREFIX}_LIBRARIES "${${acl_PREFIX}_LIBRARIES}"  CACHE INTERNAL "List of component targets" FORCE)
      set( ${acl_PREFIX}_PACKAGE_LIST "${${acl_PREFIX}_PACKAGE_LIST}"  CACHE INTERNAL
         "List of known ${acl_PREFIX} targets" FORCE)
      set( ${acl_PREFIX}_TPL_LIST "${${acl_PREFIX}_TPL_LIST}"  CACHE INTERNAL 
         "List of third party libraries known by ${acl_PREFIX}" FORCE)
      set( ${acl_PREFIX}_TPL_INCLUDE_DIRS "${${acl_PREFIX}_TPL_INCLUDE_DIRS}"  CACHE
         INTERNAL "List of include paths used by ${acl_PREFIX} to find thrid party vendor header files." 
         FORCE)
      set( ${acl_PREFIX}_TPL_LIBRARIES "${${acl_PREFIX}_TPL_LIBRARIES}"  CACHE INTERNAL
         "List of third party libraries used by ${acl_PREFIX}." FORCE)
      set( ${acl_PREFIX}_EXPORT_TARGET_PROPERTIES
         "${${acl_PREFIX}_EXPORT_TARGET_PROPERTIES}" PARENT_SCOPE)

   endif()

  
endmacro()

# ------------------------------------------------------------
# Register_scalar_test()
#
# 1. Special treatment for Roadrunner/ppe code (must ssh and then run)
# 2. Register the test
# 3. Register the pass/fail criteria.
# ------------------------------------------------------------
macro( register_scalar_test targetname runcmd command cmd_args )
   
   # Cielito needs the ./ in front of the binary name.
   if( "${MPIEXEC}" MATCHES "aprun" OR "${MPIEXEC}" MATCHES "srun" )
      set( APT_TARGET_FILE_PREFIX "./" )
   endif()
   separate_arguments( cmdargs UNIX_COMMAND ${cmd_args} )
   add_test( 
      NAME    ${targetname}
      COMMAND ${RUN_CMD} ${APT_TARGET_FILE_PREFIX}${command} ${cmdargs}
   )

   # reserve enough threads for application unit tests
   set( num_procs 1 ) # normally we only need 1 core for each scalar
                      # test.  For application unit tests, a parallel
                      # job is forked that needs more cores.
   if( addscalartest_APPLICATION_UNIT_TEST )
      if( "${cmd_args}" MATCHES "--np" AND NOT "${cmd_args}" MATCHES "scalar")
         string( REGEX REPLACE "--np ([0-9]+)" "\\1" num_procs
                      "${cmd_args}" )
         # the forked processes needs $num_proc threads.  add one for
         # the master thread, the original scalar process.
         math( EXPR num_procs  "${num_procs} + 1" )
      endif()
   endif()

   # set pass fail criteria, processors required, etc.
   set_tests_properties( ${targetname}
      PROPERTIES	
      PASS_REGULAR_EXPRESSION "${addscalartest_PASS_REGEX}"
      FAIL_REGULAR_EXPRESSION "${addscalartest_FAIL_REGEX}"
      PROCESSORS              "${num_procs}"
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

   # Labels
   # message("LABEL (ast) = ${addscalartest_LABEL}")
   if( NOT "${addscalartest_LABEL}x" STREQUAL "x" )
      set_tests_properties( ${targetname}
         PROPERTIES  LABELS "${addscalartest_LABEL}" )
   endif()
endmacro()

# ------------------------------------------------------------
# Register_parallel_test()
#
# 1. Register the test
# 2. Register the pass/fail criteria.
# ------------------------------------------------------------
macro( register_parallel_test targetname numPE command cmd_args )
   if( VERBOSE )
      message( "      Adding test: ${targetname}" )
   endif()
   if( addparalleltest_MPI_PLUS_OMP )
      add_test( 
         NAME    ${targetname}
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${numPE}
                 ${MPIEXEC_OMP_POSTFLAGS}
                 ${command}
                 ${cmdarg}
                 )
   else()
      add_test( 
         NAME    ${targetname}
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${numPE}
                 ${MPIRUN_POSTFLAGS}
                 ${command}
                 ${cmdarg}
                 )
   endif()
   set_tests_properties( ${targetname}
      PROPERTIES	
        PASS_REGULAR_EXPRESSION "${addparalleltest_PASS_REGEX}"
        FAIL_REGULAR_EXPRESSION "${addparalleltest_FAIL_REGEX}"
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
   if( addparalleltest_MPI_PLUS_OMP ) 
      math( EXPR numthreads "${numPE} * ${MPI_CORES_PER_CPU}" )
      # message("target = ${targetname}, numthreads = ${numthreads}")
      if( MPI_HYPERTHREADING )
         math( EXPR numthreads "2 * ${numthreads}" )
         # message("target = ${targetname}, numthreads = ${numthreads}")
      endif()
      # message("target = ${targetname}, numthreads = ${numthreads}")
      set_tests_properties( ${targetname}
         PROPERTIES
           PROCESSORS "${numthreads}"
           # RUN_SERIAL "ON"
           LABELS "nomemcheck" )
   else()
     if( "${addparalleltest_LABEL}x" STREQUAL "x" )
       set_tests_properties( ${targetname}
         PROPERTIES PROCESSORS "${numPE}" )
     else()
       # message("LABEL (apt) = ${addparalleltest_LABEL}")
       set_tests_properties( ${targetname}
         PROPERTIES PROCESSORS "${numPE}"
         LABELS "${addparalleltest_LABEL}" )
     endif()
   endif()
endmacro()

#----------------------------------------------------------------------#
# Special post-build options for Win32 platforms
# ------------------------------------------------------------
# copy_win32_dll_to_test_dir()
#
#----------------------------------------------------------------------#
macro( copy_win32_dll_to_test_dir )
   if( WIN32 )
      # For Win32 with shared libraries, the package dll must be
      # located in the test directory.

      # Discover all library dependencies for this unit test.
      get_target_property( link_libs Ut_${compname}_${testname}_exe LINK_LIBRARIES )
      set( old_link_libs "" )
                
      # Recurse through the library dependencies to build a list of all .dll dependencies.
      while( NOT "${old_link_libs}" STREQUAL "${link_libs}" )
         set( old_link_libs ${link_libs} )
         foreach( lib ${link_libs} )
            # $lib will either be a cmake target (e.g.: Lib_dsxx, Lib_c4) or an actual path
            # to a library (c:\lib\gsl.lib).
            if( NOT EXISTS ${lib} )
               # Must be a CMake target... find it's dependencies...
               get_target_property( link_libs2 ${lib} LINK_LIBRARIES )
               list( APPEND link_libs ${link_libs2} )
            endif()
         endforeach()
         list( REMOVE_DUPLICATES link_libs )
         foreach( lib ${link_libs} )
            if( "${lib}" MATCHES "NOTFOUND" )
               # nothing to add so remove from list
               list( REMOVE_ITEM link_libs ${lib} )
            elseif( "${lib}" MATCHES ".[lL]ib$" )
               # We have a path to a static library. Static libraries do not 
               # need to be copied.  
               list( REMOVE_ITEM link_libs ${lib} )
               # However, if there is a corresponding dll, we should add it 
               # to the list.
               string( REPLACE ".lib" ".dll" dll_lib ${lib} )
               if( ${dll_lib} MATCHES "[.]dll$" AND EXISTS ${dll_lib} )
                  list( APPEND link_libs "${dll_lib}" )
               endif()               
            endif()
         endforeach()
      endwhile()
      list( REMOVE_DUPLICATES link_libs )      
      
      # Add a post-build command to copy each dll into the test directory.
      foreach( lib ${link_libs} )
         unset( ${comp_target}_loc )
         if( EXISTS ${lib} )
            # If $lib is a full path to a library, add it to the list
            set( ${comp_target}_loc ${lib} )
            set( ${comp_target}_gnutoms NOTFOUND )
         else()
            # if $lib is a target name, obtain the file path.
            set( ${comp_target}_loc $<TARGET_FILE:${lib}> )
            get_target_property( ${comp_target}_gnutoms ${lib} GNUtoMS )
         endif()
         # Also grab the file with debug info
         string( REPLACE ".dll" ".pdb" pdb_file ${${comp_target}_loc} )

         if( "${comp_target}_loc" MATCHES "rtt" AND NOT ${comp_target}_gnutoms )
            # message("copy dbg")
            add_custom_command( TARGET Ut_${compname}_${testname}_exe 
               POST_BUILD
               COMMAND ${CMAKE_COMMAND} -E copy_if_different ${${comp_target}_loc} 
                       ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
               COMMAND ${CMAKE_COMMAND} -E copy_if_different ${pdb_file} 
                       ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
               )
         else()
            add_custom_command( TARGET Ut_${compname}_${testname}_exe 
               POST_BUILD
               COMMAND ${CMAKE_COMMAND} -E copy_if_different ${${comp_target}_loc} 
                       ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
               )
         endif()
      endforeach()
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
# add_scalar_tests( 
#    SOURCES "${test_sources}" 
#    [ DEPS    "${library_dependencies}" ]
#    [ TEST_ARGS     "arg1;arg2" ]
#    [ PASS_REGEX    "regex" ]
#    [ FAIL_REGEX    "regex" ]
#    [ RESOURCE_LOCK "lockname" ]
#    [ RUN_AFTER     "test_name" ]
#    [ APPLICATION_UNIT_TEST ]
# )
#
# Options:
#   APPLICATION_UNIT_TEST - (CI/CT only) If present, do not run the
#        test under 'aprun'.  ApplicationUnitTest based tests must be
#        run this way.  Setting this option when DRACO_C4==SCALAR will
#        reset any value provided in TEST_ARGS to be "--np scalar".
#   LINK_WITH_FORTRAN - Tell the compiler to use the Fortran compiler
#        for the final link of the test.  This is needed for Intel and
#        PGI. 
#
#----------------------------------------------------------------------#
macro( add_scalar_tests test_sources )

   # These become variables of the form ${addscalartests_SOURCES}, etc.
   parse_arguments( 
      # prefix
      addscalartest
      # list names
      "SOURCES;DEPS;TEST_ARGS;PASS_REGEX;FAIL_REGEX;RESOURCE_LOCK;RUN_AFTER;LABEL"
      # option names
      "APPLICATION_UNIT_TEST;LINK_WITH_FORTRAN;NONE"
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
   elseif( "${MPIEXEC}" MATCHES "srun" )
      set( RUN_CMD ${MPIEXEC} -n 1 )
   else()
      unset( RUN_CMD )
   endif()

   # Special cases for tests that use the ApplicationUnitTest
   # framework (see c4/ApplicationUnitTest.hh).
   if( addscalartest_APPLICATION_UNIT_TEST )
      # This is a special case for catamount (CI/CT).  For appliation
      # unit tests, the main test runs on the 'login' node (1 rank
      # only) and the real test is run under 'aprun'.  So we do not
      # prefix the test command with 'aprun'.
      if( "${MPIEXEC}" MATCHES "aprun" )
         unset( RUN_CMD )
      endif()
      
      # If this is an ApplicationUnitTest based test then the
      # TEST_ARGS will look like "--np 1;--np 2;--np 4".  For the case
      # where DRACO_C4 = SCALAR, we will automatically demote these
      # arguments to "--np scalar."
      if( "${DRACO_C4}" MATCHES "SCALAR" )
         set( addscalartest_TEST_ARGS "--np scalar" )
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

   # Loop over each test source files:
   # 1. Compile the executable
   # 2. Register the unit test

   # Generate the executable
   # ------------------------------------------------------------
   foreach( file ${addscalartest_SOURCES} )
   
      if( "${file}" MATCHES "tstParallelUnitTest" )
         message("RUN_CMD = ${RUN_CMD}")
      endif()

      get_filename_component( testname ${file} NAME_WE )
      add_executable( Ut_${compname}_${testname}_exe ${file} )
      set_target_properties( Ut_${compname}_${testname}_exe 
         PROPERTIES 
           OUTPUT_NAME ${testname} 
           VS_KEYWORD  ${testname}
           FOLDER      ${compname}_test
         )
      # Do we need to use the Fortran compiler as the linker?
      if( addscalartest_LINK_WITH_FORTRAN )
         set_target_properties( Ut_${compname}_${testname}_exe 
            PROPERTIES LINKER_LANGUAGE Fortran )
      endif()
      target_link_libraries( 
         Ut_${compname}_${testname}_exe 
         ${test_lib_target_name}
         ${addscalartest_DEPS}
         )         
         
      # Special post-build options for Win32 platforms
      # ------------------------------------------------------------
      copy_win32_dll_to_test_dir()
    
   endforeach()
         
   # Register the unit test
   # ------------------------------------------------------------
   foreach( file ${addscalartest_SOURCES} )
      get_filename_component( testname ${file} NAME_WE )

      if( "${addscalartest_TEST_ARGS}none" STREQUAL "none" )
         register_scalar_test( ${compname}_${testname} 
            "${RUN_CMD}" ${testname} "" )
       else()
          set( iarg "0" )
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
      "SOURCES;PE_LIST;DEPS;TEST_ARGS;PASS_REGEX;FAIL_REGEX;RESOURCE_LOCK;RUN_AFTER;MPIFLAGS;LABEL"
      # option names
      "MPI_PLUS_OMP;LINK_WITH_FORTRAN"
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
      list( APPEND addparalleltest_FAIL_REGEX ".*ERROR:.*" )
   endif()

   # Format resource lock command
   if( NOT "${addparalleltest_RESOURCE_LOCK}none" STREQUAL "none" )
      set( addparalleltest_RESOURCE_LOCK 
         "RESOURCE_LOCK ${addparalleltest_RESOURCE_LOCK}")
   endif()

   # What is the component name? Use this to give a target name to the test.
   string( REPLACE "_test" "" compname ${PROJECT_NAME} )

   # Override MPI Flags upon user request
   if ( NOT DEFINED addparalleltest_MPIFLAGS )
       set( MPIRUN_POSTFLAGS ${MPIEXEC_POSTFLAGS} )
   else()
       set( MPIRUN_POSTFLAGS "${addparalleltest_MPIFLAGS}" )
   endif()
  

   # Loop over each test source files:
   # 1. Compile the executable
   # 2. Link against dependencies (libraries)

   foreach( file ${addparalleltest_SOURCES} )
      get_filename_component( testname ${file} NAME_WE )
      if( VERBOSE )
         message( "   add_executable( Ut_${compname}_${testname}_exe ${file} )")
      endif()
      add_executable( Ut_${compname}_${testname}_exe ${file} )
      set_target_properties( 
         Ut_${compname}_${testname}_exe 
         PROPERTIES 
           OUTPUT_NAME ${testname} 
           VS_KEYWORD  ${testname}
           FOLDER      ${compname}_test
         )
      if( addparalleltest_MPI_PLUS_OMP )
         if( ${CMAKE_GENERATOR} MATCHES Xcode )
            set_target_properties( Ut_${compname}_${testname}_exe 
               PROPERTIES XCODE_ATTRIBUTE_ENABLE_OPENMP_SUPPORT YES )
         endif()
      endif()
      # Do we need to use the Fortran compiler as the linker?
      if( addparalleltest_LINK_WITH_FORTRAN )
           set_target_properties( Ut_${compname}_${testname}_exe 
               PROPERTIES LINKER_LANGUAGE Fortran )
      endif()

      target_link_libraries( 
         Ut_${compname}_${testname}_exe 
         ${test_lib_target_name} 
         ${addparalleltest_DEPS}
         )

      # Special post-build options for Win32 platforms
      # ------------------------------------------------------------
      copy_win32_dll_to_test_dir()
         
   endforeach()

   # 3. Register the unit test
   # 4. Register the pass/fail criteria.
   if( ${DRACO_C4} MATCHES "MPI" )
      foreach( file ${addparalleltest_SOURCES} )

         get_filename_component( testname ${file} NAME_WE )
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
# provide_aux_files
# 
# Call this macro from a package CMakeLists.txt to instruct the build 
# system that some files should be copied from the source directory
# into the build directory.
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
      set( outfile ${PROJECT_BINARY_DIR}/${srcfilenameonly} )
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
    set( folder_name ${compname}_test )
    set_target_properties( Ut_${compname}_install_inputs_${Ut_${compname}_install_inputs_iarg} 
      PROPERTIES FOLDER ${folder_name}
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
# PROCESS_AUTODOC_PAGES - Run configure_file(...) for all .dcc.in
# files found in the autodoc directory.  Destination will be the
# autodoc directory in the component binary directory.  The
# CMakeLists.txt in the draco/autodoc directory knows how to find
# these files.
#
# This allows CMAKE variables to be inserted into the .dcc files
# (e.g.: @DRACO_VERSION@)
# 
# E.g.: process_autodoc_pages()
#----------------------------------------------------------------------#
macro( process_autodoc_pages )
   file( GLOB autodoc_in autodoc/*.in )
   foreach( file ${autodoc_in} )
      get_filename_component( dest_file ${file} NAME_WE )
      configure_file( ${file} ${PROJECT_BINARY_DIR}/autodoc/${dest_file}.dcc @ONLY )
   endforeach()
endmacro()

#----------------------------------------------------------------------#
# End
#----------------------------------------------------------------------#
