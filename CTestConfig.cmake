## This file should be placed in the root directory of your project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.
## # The following are required to uses Dart and the Cdash dashboard
##   ENABLE_TESTING()
##   INCLUDE(CTest)
set( CTEST_PROJECT_NAME "Draco" )
set( CTEST_NIGHTLY_START_TIME "21:00:00 MDT" )

set( CTEST_DROP_METHOD "https")
set( CTEST_DROP_SITE "ccstest.lanl.gov" )
set( CTEST_DROP_LOCATION  "/cdash/submit.php?project=${CTEST_PROJECT_NAME}" )
set( CTEST_DROP_SITE_CDASH TRUE )
set( CTEST_CURL_OPTIONS CURLOPT_SSL_VERIFYPEER_OFF CURLOPT_SSL_VERIFYHOST_OFF )
