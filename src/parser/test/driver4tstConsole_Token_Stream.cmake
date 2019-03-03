#-----------------------------*-cmake-*----------------------------------------#
# file   parser/test/driver4tstConsole_Token_Stream.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   Friday, Jul 10, 2015, 14:16 pm
# brief  This is a CTest script that is used to test parser/Ipcress_Interpreter
# note   Copyright (C) 2016, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id: CMakeLists.txt 6721 2012-08-30 20:38:59Z gaber $
#------------------------------------------------------------------------------#

# See examples at config/ApplicationUnitTest.cmake and diagnostics/test/tDracoInfo.cmake

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
## Check the output
##---------------------------------------------------------------------------##

message("Checking the generated output file...
")
if(WIN32)
  set( exe_suffix ".exe")
endif()

# This string should be found:

string(FIND "${testout}" "tstConsole_Token_Stream${exe_suffix} Test: PASSED." string_pos)
if( ${string_pos} GREATER 0 )
  PASSMSG( "tstConsole_Token_Stream ran successfully." )
else()
  FAILMSG( "tstConsole_Token_Stream did not run successfully." )
endif()

##---------------------------------------------------------------------------##
## Final report
##---------------------------------------------------------------------------##
aut_report()

##---------------------------------------------------------------------------##
## End
##---------------------------------------------------------------------------##
