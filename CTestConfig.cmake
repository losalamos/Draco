#-----------------------------*-cmake-*----------------------------------------#
# file   CTestConfig.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 Nov 22
# brief  Link to CDash build/test results dashboard
# note   Copyright © 2010 Los Alamos National Security
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

## This file should be placed in the root directory of your project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.
## The following are required to uses Dart and the Cdash dashboard
##   ENABLE_TESTING()
##   INCLUDE(CTest)
## http://www.cmake.org/Wiki/CMake_Testing_With_CTest

if( NOT CTEST_PROJECT_NAME )
  set( CTEST_PROJECT_NAME "Draco")
endif()
set( CTEST_NIGHTLY_START_TIME "00:00:00 MST")

set( CTEST_DROP_METHOD "http")
set( CTEST_DROP_SITE "coder.lanl.gov")
set( CTEST_DROP_LOCATION 
   "/cdash/submit.php?project=${CTEST_PROJECT_NAME}" )
set( CTEST_DROP_SITE_CDASH TRUE )
set( CTEST_CURL_OPTIONS CURLOPT_SSL_VERIFYPEER_OFF )

# Options:
# set(UPDATE_TYPE "true")
# set(MEMORYCHECK_SUPPRESSIONS_FILE ${CMAKE_SOURCE_DIR}/tests/valgrind-csync.supp)
# set( CTEST_TRIGGER_SITE 
#      "http://${DROP_SITE}/cgi-bin/Submit-CMake-TestingResults.pl" )
