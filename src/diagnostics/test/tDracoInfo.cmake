#-----------------------------*-cmake-*----------------------------------------#
# file   diagnostics/test/tDracoInfo.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   Monday, Nov 19, 2012, 17:02 pm
# brief  This is a CTest script that is used to test bin/draco_info.
# note   Copyright (C) 2012-2015, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id: CMakeLists.txt 6721 2012-08-30 20:38:59Z gaber $
#------------------------------------------------------------------------------#
#

# Some useful macros
get_filename_component( draco_config_dir
   ${CMAKE_CURRENT_LIST_DIR}/../../../config ABSOLUTE )
set( CMAKE_MODULE_PATH ${draco_config_dir} )
include( ApplicationUnitTest )

# Setup and Sanity check provides:
#   APP       - cleaned up path to executable
#   OUTFILE   - Output filename
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

message("${RUN_CMD} ${APP} > ${OUTFILE}")
separate_arguments(RUN_CMD)
execute_process(
  COMMAND ${RUN_CMD} ${APP}
  WORKING_DIRECTORY ${WORKDIR}
  RESULT_VARIABLE testres
  OUTPUT_VARIABLE testout
  ERROR_VARIABLE  testerror
  )

##---------------------------------------------------------------------------##
# Ensure there are no errors
if( NOT "${testres}" STREQUAL "0" )
   message( FATAL_ERROR "Test FAILED:
     error message = ${testerror}")
endif()

##---------------------------------------------------------------------------##
## Echo the output to stdout and to an output file for parsing.
if( VERBOSE )
   message("${testout}")
endif()
file( WRITE ${OUTFILE} ${testout} )

##---------------------------------------------------------------------------##
## Analyize the output directly.
##---------------------------------------------------------------------------##

string( REGEX REPLACE "\n" ";" testout ${testout} )
set( foundcopyright FALSE )
set( foundsystemtype FALSE )
foreach( line ${testout} )
   if( ${line} MATCHES "Copyright [(]C[)]" )
      set( foundcopyright TRUE )
   endif()
   if( ${line} MATCHES "System type" )
      set( foundsystemtype TRUE )
   endif()

   #    string( REGEX REPLACE ".*= ([0-9.]+).*" "\\1" value ${line} )
   #    set( refvalue "6411.71" )
   #    if( ${refvalue} EQUAL ${value} )
   #       PASSMSG( "Specific Ion Internal Energy matches expected value.")
   #    else()
   #       FAILMSG( "Specific Ion Internal Energy does not match expected value.")
   #    endif()
   # endif()
endforeach()

if( foundcopyright )
   PASSMSG( "Found copyright date")
else()
   FAILMSG( "Did not find copyright date")
endif()
if( foundsystemtype )
   PASSMSG( "Found system type id")
else()
   FAILMSG( "Did not find system type id")
endif()

##---------------------------------------------------------------------------##
## End
##---------------------------------------------------------------------------##
