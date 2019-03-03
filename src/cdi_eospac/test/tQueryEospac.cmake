#-----------------------------*-cmake-*----------------------------------------#
# file   cdi_eospac/test/tQueryEospac.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   Friday, Jul 10, 2015, 14:16 pm
# brief  This is a CTest script that is used to test cdi_eospac/QueryEospac
# note   Copyright (C) 2016, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

# Setup the CMake based ApplicationUnitTest environment
set( CMAKE_MODULE_PATH ${DRACO_CONFIG_DIR} )
include( ApplicationUnitTest )

# Setup and Sanity check provides:
aut_setup()

##---------------------------------------------------------------------------##
# Run the application and capture the output.
# Variables available for inspection:
#   ${testres} contains the return code
#   ${testout} contains stdout
#   ${testerror} contains stderr
aut_runTests()

##---------------------------------------------------------------------------##
## Case 1: Analyze the output vs. a gold file
##---------------------------------------------------------------------------##
if( GOLDFILE )
  aut_numdiff()

   ## Analyze the output directly.
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

##---------------------------------------------------------------------------##
## Check output for --version and --help versions.
##---------------------------------------------------------------------------##
else()

  if( ${ARGVALUE} STREQUAL "--version" )
    string(FIND "${testout}" "QueryEospac: version" POS1)
    string(FIND "${testout}" "QueryEospac.exe: version" POS2)
    if( ${POS1} GREATER 0 OR ${POS2} GREATER 0 )
      PASSMSG( "Version tag found in the output.")
    else()
      FAILMSG( "Version tag NOT found in the output.")
    endif()

  elseif(${ARGVALUE} STREQUAL "--help" )
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
## Final report
##---------------------------------------------------------------------------##
aut_report()

##---------------------------------------------------------------------------##
## End
##---------------------------------------------------------------------------##
