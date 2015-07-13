#-----------------------------*-cmake-*----------------------------------------#
# file   config/ApplicationUnitTest.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   Monday, Nov 19, 2012, 16:21 pm
# brief  Provide macros that aid in creating unit tests that run
#        interactive user codes (i.e.: run a binary that reads an
#        input file and diff the resulting output file).
# note   Copyright (C) 2012-2015, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id: CMakeLists.txt 6732 2012-09-05 22:28:18Z kellyt $
#------------------------------------------------------------------------------#

# Reference: https://rtt.lanl.gov/redmine/projects/draco/wiki/CMake-based_ApplicationUnitTest
# Example: draco/src/diagnostics/test/tDracoInfo.cmake (and associated CMakeLists.txt).

include( parse_arguments )

##---------------------------------------------------------------------------##
## Check values for $APP
##
## Requires:
##  APP        - name of executable that will be run
##
## Returns:
##   APP       - cleaned up path to executable
##   STDINFILE - cleaned up path to intput file.
##   BINDIR    - Directory location of binary file
##   PROJECT_BINARY_DIR - Parent directory of BINDIR
##   GOLDFILE  - cleaned up path to gold standard file.
##   OUTFILE   - Output filename derived from GOLDFILE.
##   ERRFILE   - Output(error) filename derived from GOLDFILE.
##
## if VERBOSE is set, also echo input values.
##   APP       - path name for executable.
##   STDINFILE - optional input file
##   GOLDFILE  - optional gold standard file.
##---------------------------------------------------------------------------##
macro( aut_setup)

   if( VERBOSE )
      message("Running tQueryEospac.cmake with the following parameters:")
      message("   APP       = ${APP}")
      if( STDINFILE )
         message("   STDINFILE = ${STDINFILE}")
      endif()
      if( GOLDFILE )
         message("   GOLDFILE = ${GOLDFILE}" )
      endif()
   endif()

   # Setup and sanity check

   if( "${APP}x" STREQUAL "x" )
      message( FATAL_ERROR "You must provide a value for APP." )
   endif()

   get_filename_component( APP ${APP} ABSOLUTE )
   if( STDINFILE )
      get_filename_component( STDINFILE ${STDINFILE} ABSOLUTE )
   endif()
   get_filename_component( BINDIR ${APP} PATH )
   get_filename_component( PROJECT_BINARY_DIR ${BINDIR} PATH )
   if( GOLDFILE )
      get_filename_component( OUTFILE ${GOLDFILE} NAME_WE )
   else()
      get_filename_component( OUTFILE ${APP} NAME_WE )
   endif()
   set( OUTFILE "${CMAKE_CURRENT_BINARY_DIR}/${OUTFILE}.out")
   set( ERRFILE "${CMAKE_CURRENT_BINARY_DIR}/${OUTFILE}.err")

   message("Testing ${APP}")

   if( NOT EXISTS "${APP}" )
      message( FATAL_ERROR "Cannot find ${APP}")
   endif()

   set( numpasses 0 )
   set( numfails  0 )

endmacro()

##---------------------------------------------------------------------------##
## PASSMSG/FAILMSG
##---------------------------------------------------------------------------##

macro(PASSMSG msg)
    math( EXPR numpasses "${numpasses} + 1" )
    message( "Test Passes: ${msg}")
endmacro()

macro(FAILMSG msg)
    math( EXPR numfails "${numfails} + 1" )
    message( "Test Fails: ${msg}")
endmacro()

macro(ITFAILS)
    math( EXPR numfails "${numfails} + 1" )
endmacro()

##---------------------------------------------------------------------------##
## REGISTRATION
##---------------------------------------------------------------------------##

macro( add_app_unit_test )

  # These become variables of the form ${addscalartests_SOURCES}, etc.
  parse_arguments(
    # prefix
    aut
    # list names
    "DRIVER;APP;WORKDIR;TEST_ARGS;PASS_REGEX;FAIL_REGEX;RESOURCE_LOCK;RUN_AFTER;LABELS"
    # option names
    "NONE"
    ${ARGV}
    )

  #
  # Check required intputs and set default vaules:
  #
  if( NOT EXISTS ${aut_DRIVER} )
    message( FATAL_ERROR "Could not find the cmake driver script = ${aut_DRIVER}." )
  endif()
  if( NOT DEFINED aut_APP )
    message( FATAL_ERROR "You must provide a value for APP." )
  endif()
  if( NOT DEFINED aut_WORKDIR )
    set( aut_WORKDIR ${CMAKE_CURRENT_BINARY_DIR} )
  endif()
  if( NOT DEFINED aut_PASS_REGEX )
    set( aut_PASS_REGEX "PASSED" )
  endif()
  if( NOT DEFINED aut_FAIL_REGEX )
    set( aut_FAIL_REGEX "FAILED;fails" )
  endif()
  if( DEFINED aut_LABELS )
    set( LABEL "LABELS ${aut_LABELS}" )
  endif()

  # Load some information from the build environment:
  unset( RUN_CMD )
  if( "${C4_MPICMD}" MATCHES "aprun" )
    set( RUN_CMD "aprun -n 1" )
  elseif( HAVE_MIC )
    set( RUN_CMD "ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null $ENV{HOSTNAME}-mic0 ${Draco_BINARY_DIR}/config/run_test_on_mic.sh ${WORKDIR}" )
  endif()

  # Create a name for the test
  get_filename_component( drivername ${aut_DRIVER} NAME_WE )
  get_filename_component( package_name ${aut_DRIVER} PATH )
  get_filename_component( package_name ${package_name} PATH )
  get_filename_component( package_name ${package_name} NAME )
  set( ctestname_base ${package_name}_${drivername} )

message("

add_test(
   NAME ${ctestname_base}
   COMMAND ${CMAKE_COMMAND}
     -D APP              = ${aut_APP}
     -D WORKDIR          = ${aut_WORKDIR}
     -D TESTNAME         = ${ctestname_base}
     -D DRACO_CONFIG_DIR = ${Draco_SOURCE_DIR}/config
     -D DRACO_INFO       = $<TARGET_FILE_DIR:Exe_draco_info>/$<TARGET_FILE_NAME:Exe_draco_info>
     -D RUN_CMD          = ${RUN_CMD}
     -P ${aut_DRIVER}
   )
set_tests_properties( diagnostics_draco_info
   PROPERTIES
     PASS_REGULAR_EXPRESSION ${aut_PASS_REGEX}
     FAIL_REGULAR_EXPRESSION ${aut_FAIL_REGEX}
     ${LABELS}
   )

")

# Register the test...

add_test(
   NAME ${ctestname_base}
   COMMAND ${CMAKE_COMMAND}
     -D APP=${aut_APP}
     -D WORKDIR=${aut_WORKDIR}
     -D TESTNAME=${ctestname_base}
     -D DRACO_CONFIG_DIR=${Draco_SOURCE_DIR}/config
     -D DRACO_INFO=$<TARGET_FILE_DIR:Exe_draco_info>/$<TARGET_FILE_NAME:Exe_draco_info>
     -D RUN_CMD=${RUN_CMD}
     -P ${aut_DRIVER}
   )
set_tests_properties( ${ctestname_base}
   PROPERTIES
     PASS_REGULAR_EXPRESSION ${aut_PASS_REGEX}
     FAIL_REGULAR_EXPRESSION ${aut_FAIL_REGEX}
     ${LABELS}
   )

endmacro()

##---------------------------------------------------------------------------##
## Run the tests
##---------------------------------------------------------------------------##
macro( aut_runTests )

  # Run the application and capture the output.
  message("
=============================================
=== ${TESTNAME}
=============================================
")
# === CMake driven ApplicationUnitTest: ${TESTNAME}

  if( DEFINED RUN_CMD )
    message(">>> Running: ${RUN_CMD} ${APP}
>>>      > ${OUTFILE}
")
  else()
    message(">>> Running: ${APP}
>>>      > ${OUTFILE}
")
  endif()

  separate_arguments(RUN_CMD)

  # Print version information
  if( EXISTS ${DRACO_INFO} )
    execute_process(
      COMMAND ${RUN_CMD} ${APP} --version
      WORKING_DIRECTORY ${WORKDIR}
      RESULT_VARIABLE testres
      OUTPUT_VARIABLE testout
      ERROR_VARIABLE  testerror
      )
    if( NOT ${testres} STREQUAL "0" )
      FAILMSG("Unable to run 'draco_info --version'")
    else()
      message("${testout}")
    endif()
  endif()

  # Run the application capturing all output.
  execute_process(
    COMMAND ${RUN_CMD} ${APP}
    WORKING_DIRECTORY ${WORKDIR}
    RESULT_VARIABLE testres
    OUTPUT_VARIABLE testout
    ERROR_VARIABLE  testerror
    )

  # Capture all the output to log files:
  file( WRITE ${OUTFILE} ${testout} )
  file( WRITE ${ERRFILE} ${testerr} )

  # Ensure there are no errors
  if( NOT "${testres}" STREQUAL "0" )
    message( FATAL_ERROR "Test FAILED:
     error message = ${testerror}")
  else()
    message("${testout}")
  endif()

endmacro()


##---------------------------------------------------------------------------##
## Run the tests
##---------------------------------------------------------------------------##
macro( aut_report )

  message("
*********************************************")
  if( ${numpasses} GREATER 0 AND ${numfails} STREQUAL 0 )
    message("**** ${TESTNAME}: PASSED.")
  else()
    message("**** ${TESTNAME}: FAILED.")
  endif()
  message("*********************************************
")

endmacro()

##---------------------------------------------------------------------------##
## end
##---------------------------------------------------------------------------##
