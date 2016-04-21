#-----------------------------*-cmake-*----------------------------------------#
# file   cdi_ipcress/test/tIpcress_Interpreter.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   Friday, Jul 10, 2015, 14:16 pm
# brief  This is a CTest script that is used to test cdi_ipcress/Ipcress_Interpreter
# note   Copyright (C) 2016, Los Alamos National Security, LLC.
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

# These strings should be found:

string(FIND "${testout}" "This opacity file has 2 materials:" string_pos)
if( ${string_pos} GREATER 0 )
  PASSMSG( "Found 2 materials." )
else()
  FAILMSG( "Did not find 2 materials." )
endif()

string(FIND "${testout}" "Material 1 has ID number 10001" string_pos)
if( ${string_pos} GREATER 0 )
  PASSMSG( "Found material ID 10001." )
else()
  FAILMSG( "Did not find material ID 10001." )
endif()

string(FIND "${testout}" "Frequency grid" string_pos)
if( ${string_pos} GREATER 0 )
  PASSMSG( "Found Frequency grid." )
else()
  FAILMSG( "Did not find Frequency grid." )
endif()
message(" ")

# Diff the output vs a gold file.
aut_numdiff()

##---------------------------------------------------------------------------##
## Final report
##---------------------------------------------------------------------------##
aut_report()

##---------------------------------------------------------------------------##
## End
##---------------------------------------------------------------------------##
