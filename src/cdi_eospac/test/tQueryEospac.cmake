#-----------------------------*-cmake-*----------------------------------------#
# This is a CTest script that is used to test bin/QueryEospac.
#
# Calling example
# cmake \
# -DAPP=QueryEospac \
# [-DCMDARGS=--version|--help] \
# [-DSTDINFILE=QueryEospac.input] \
# [-DGOLDFILE=QueryEospac.gold] \
# -DDraco_BINARY_DIR=${Draco_BINARY_DIR} \
# -P tQueryEospac.cmake

# Some useful macros
get_filename_component( draco_config_dir
   ${CMAKE_CURRENT_LIST_DIR}/../../../config ABSOLUTE )
set( CMAKE_MODULE_PATH ${draco_config_dir} )
include( ApplicationUnitTest )

# Setup and Sanity check provides:
#   APP       - cleaned up path to executable
#   STDINFILE - cleaned up path to intput file.
#   BINDIR    - Directory location of binary file
#   PROJECT_BINARY_DIR - Parent directory of BINDIR
#   GOLDFILE  - cleaned up path to gold standard file.
#   OUTFILE   - Output filename derived from GOLDFILE.
aut_setup()

##---------------------------------------------------------------------------##
# Run the application and capture the output.
message("Running tests...")

unset( RUN_CMD )
file( STRINGS ${Draco_BINARY_DIR}/src/c4/c4/config.h C4_MPICMD REGEX C4_MPICMD )
if( "${C4_MPICMD}" MATCHES "aprun" )
  set( RUN_CMD "aprun -n 1" )
elseif( HAVE_MIC )
  set( RUN_CMD "ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null $ENV{HOSTNAME}-mic0 ${Draco_BINARY_DIR}/config/run_test_on_mic.sh ${WORKDIR}" )
endif()

##---------------------------------------------------------------------------##
## Case 1: Read commands from an input file.
if( DEFINED GOLDFILE )
  message( "${RUN_CMD} ${APP}
   < ${STDINFILE}
   > ${OUTFILE}")
  separate_arguments(RUN_CMD)
  execute_process(
    COMMAND ${RUN_CMD} ${APP}
    INPUT_FILE ${STDINFILE}
    RESULT_VARIABLE testres
    OUTPUT_VARIABLE testout
    ERROR_VARIABLE  testerror
    )

##---------------------------------------------------------------------------##
## Case 2: test command line option --version or --help
elseif( DEFINED CMDARGS )
  message( "${RUN_CMD} ${APP} ${CMDARGS}")
  separate_arguments(RUN_CMD)
  execute_process(
    COMMAND ${RUN_CMD} ${APP} ${CMDARGS}
    RESULT_VARIABLE testres
    OUTPUT_VARIABLE testout
    ERROR_VARIABLE  testerror
    )
  set( OUTFILE "QueryEospac${CMDARGS}.out" )

else()
  message( FATAL_ERROR "expected CMDARGS or GOLDFILE to be defined.")
endif()

##---------------------------------------------------------------------------##
# Ensure there are no errors
if( NOT "${testerror}x" STREQUAL "x" OR NOT "${testres}" STREQUAL "0" )
   message( FATAL_ERROR "Test FAILED:
     mesage = ${testerror}")
endif()

##---------------------------------------------------------------------------##
## Echo the output to stdout and to an output file for parsing.
if( VERBOSE )
   message("${testout}")
endif()
file( WRITE ${OUTFILE} ${testout} )

##---------------------------------------------------------------------------##
## Analyze the output vs. a gold file
if( GOLDFILE )
   ## Use numdiff to compare output
   find_program( exenumdiff numdiff )
   if( NOT EXISTS ${exenumdiff} )
      message( FATAL_ERROR "Numdiff not found in PATH")
   endif()
   message("${exenumdiff}
   ${OUTFILE}
   ${GOLDFILE}")
   execute_process(
      COMMAND ${exenumdiff} ${OUTFILE} ${GOLDFILE}
      RESULT_VARIABLE numdiffres
      OUTPUT_VARIABLE numdiffout
      ERROR_VARIABLE numdifferror
      )
   if( ${numdiffres} STREQUAL 0 )
      PASSMSG("gold matches out.")
   else()
      FAILMSG("gold does not match out.
numdiff output = ${numdiffout}" )
   endif()

   ## Analyize the output directly.
   string( REGEX REPLACE "\n" ";" testout ${testout} )
   foreach( line ${testout} )
     if( ${line} MATCHES "Specific Ion Internal Energy" )
       string( REGEX REPLACE ".*= ([0-9.]+).*" "\\1" value ${line} )
       set( refvalue "6411.71" )
       if( ${refvalue} EQUAL ${value} )
         PASSMSG( "Specific Ion Internal Energy matches expected value.")
       else()
         FAILMSG( "Specific Ion Internal Energy does not match expected value.")
       endif()
     endif()
   endforeach()

else()

  # CMDARGS versions...

  if( ${CMDARGS} STREQUAL "--version" )
    string(FIND "${testout}" "QueryEospac: version" POS)
    if(  ${POS} GREATER 0 )
      PASSMSG( "Version tag found in the output.")
    else()
      FAILMSG( "Version tag NOT found in the output.")
    endif()

  elseif(${CMDARGS} STREQUAL "--help" )
    string(FIND "${testout}"
      "Follow the prompts to print equation-of-state data to the screen."
      POS REVERSE)
    if( ${POS} GREATER 0 )
      PASSMSG( "Help prompt was found in the output." )
    else()
      FAILMSG( "Help prompt was NOT found in the output." )
    endif()
  endif()

endif()

##---------------------------------------------------------------------------##
## End
##---------------------------------------------------------------------------##
