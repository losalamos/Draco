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
set( founda FALSE )
set( foundb FALSE )
set( foundc FALSE )
set( foundinvalid FALSE )
foreach( line ${testout} )

   if( ${line} MATCHES "aflag = 1" )
      set( founda TRUE )
   endif()
   if( ${line} MATCHES "bflag = 1" )
      set( foundb TRUE )
   endif()
   if( ${line} MATCHES "cvalue = fish" )
      set( foundcfish TRUE )
   endif()
   if( ${line} MATCHES "invalid option" )
      # foo: invalid option -- 'X'
      # Usage: make [options] [target] ...
      # Options:
      #   -b, -m
      set( foundinvalid TRUE )
   endif()

endforeach()

# There are 3 versions of this test

if( ARGVALUE )

  if( "${ARGVALUE}" STREQUAL "-a" OR "${ARGVALUE}" STREQUAL "--add" )
    if( founda )
      PASSMSG( "Found a")
    else()
      FAILMSG( "Did not find a")
    endif()
  endif()
  if( ${ARGVALUE} STREQUAL "-b" )
    if( foundb )
      PASSMSG( "Found b")
    else()
      FAILMSG( "Did not find b")
    endif()
  endif()
  # This compact form is not yet supported.
  # if( ${ARGVALUE} STREQUAL "-ab" )
  #   if( foundb AND founda )
  #     PASSMSG( "Found a and b")
  #   else()
  #     FAILMSG( "Did not find a or b")
  #   endif()
  # endif()

  # KT this fails right now (need to update add_app_unit_test to allow
  # expected error return code.)  Should this print the help message?
  # if( ${ARGVALUE} STREQUAL "--badarg" )
  #   if( foundinvalid )
  #     PASSMSG( "Invalid option reported.")
  #   else()
  #     FAILMSG( "Failed to report --badarg as an invalid option.")
  #   endif()
  # endif()

  if( ${ARGVALUE} MATCHES "-c" )
    if( foundcfish )
      PASSMSG( "Found c=fish")
    else()
      FAILMSG( "Did not find c=fish")
    endif()
  endif()

else ()  # no arguments

  if( founda )
    FAILMSG( "Found a")
  else()
    PASSMSG( "Did not find a")
  endif()
   if( foundb )
    FAILMSG( "Found b")
  else()
    PASSMSG( "Did not find b")
  endif()
  if( foundc )
    FAILMSG( "Found c")
  else()
    PASSMSG( "Did not find c")
  endif()

endif()


##---------------------------------------------------------------------------##
## Final report
##---------------------------------------------------------------------------##
aut_report()

##---------------------------------------------------------------------------##
## End
##---------------------------------------------------------------------------##
