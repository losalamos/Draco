#-----------------------------*-cmake-*----------------------------------------#
# file   diagnostics/test/tDracoInfo.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   Monday, Nov 19, 2012, 17:02 pm
# brief  This is a CTest script that is used to test bin/draco_info.
# note   Copyright (C) 2016, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id: CMakeLists.txt 6721 2012-08-30 20:38:59Z gaber $
#------------------------------------------------------------------------------#

# Use config/ApplicationUnitTest.cmake test registration:
#
# include( ApplicationUnitTest )
# add_app_unit_test(
#   DRIVER ${CMAKE_CURRENT_SOURCE_DIR}/tDracoInfo.cmake
#   APP    $<TARGET_FILE_DIR:Exe_draco_info>/$<TARGET_FILE_NAME:Exe_draco_info>
#   LABELS nomemcheck )

# The above will generate a test with data similar to this:
#
# add_test(
#    NAME diagnostics_tDracoInfo
#    COMMAND /yellow/usr/projects/draco/vendors/cmake-3.2.2-Linux-x86_64/bin/cmake
#      -D APP              = $<TARGET_FILE_DIR:Exe_draco_info>/$<TARGET_FILE_NAME:Exe_draco_info>
#      -D WORKDIR          = /users/kellyt/build/ml/intel-mpid/d/src/diagnostics/test
#      -D TESTNAME         = diagnostics_tDracoInfo
#      -D DRACO_CONFIG_DIR = /users/kellyt/draco/config
#      -D DRACO_INFO       = /users/kellyt/build/ml/intel-mpid/d/src/diagnostics/draco_info
#      -D RUN_CMD          =
#      -P /users/kellyt/draco/src/diagnostics/test/tDracoInfo.cmake
#    )
# set_tests_properties( diagnostics_draco_info
#    PROPERTIES
#      PASS_REGULAR_EXPRESSION Passes
#      FAIL_REGULAR_EXPRESSION Fails
#      LABELS nomemcheck
#    )

# Variables defined above can be used in this script.

#------------------------------------------------------------------------------#
# Setup the CMake based ApplicationUnitTest environment
set( CMAKE_MODULE_PATH ${DRACO_CONFIG_DIR} )
include( ApplicationUnitTest )

# Setup and Sanity check provides:
#   APP       - cleaned up path to executable
#   OUTFILE   - Output filename
aut_setup()

##---------------------------------------------------------------------------##
# Run the application and capture the output.
# Variables available for inspection:
#   ${testres} contains the return code
#   ${testout} contains stdout
#   ${testerror} contains stderr
aut_runTests()

##---------------------------------------------------------------------------##
## Examine the output to determine if the test passed
##---------------------------------------------------------------------------##

string( REGEX REPLACE "\n" ";" testout ${testout} )
set( foundcopyright FALSE )
set( foundsystemtype FALSE )
set( foundbuilddate FALSE )
foreach( line ${testout} )
   if( ${line} MATCHES "Copyright [(]C[)]" )
      set( foundcopyright TRUE )
   endif()
   if( ${line} MATCHES "System type" )
      set( foundsystemtype TRUE )
   endif()
   if( ${line} MATCHES "build date" )
      set( foundbuilddate TRUE )
   endif()

endforeach()

# There are 3 versions of this test

if( ARGVALUE )

  if( ${ARGVALUE} STREQUAL "--version" )
    if( foundbuilddate )
      PASSMSG( "Found build date")
    else()
      FAILMSG( "Did not find build date")
    endif()
    if( foundcopyright )
      FAILMSG( "Found copyright date")
    else()
      PASSMSG( "Did not find copyright date")
    endif()

  elseif( ${ARGVALUE} STREQUAL "--brief" )
    if( foundbuilddate )
      PASSMSG( "Found build date")
    else()
      FAILMSG( "Did not find build date")
    endif()
    if( foundcopyright )
      PASSMSG( "Found copyright date")
    else()
      FAILMSG( "Did not find copyright date")
    endif()
    if( foundsystemtype )
      FAILMSG( "Found system type id")
    else()
      PASSMSG( "Did not find system type id")
    endif()

  endif()

else ()  # no arguments

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
endif()


##---------------------------------------------------------------------------##
## Final report
##---------------------------------------------------------------------------##
aut_report()

##---------------------------------------------------------------------------##
## End
##---------------------------------------------------------------------------##
