#
# Linux 64-bit
# Ref: http://www.cmake.org/Wiki/CMake_Scripting_Of_CTest

cmake_minimum_required(VERSION 3.0.0)

# Use:
# - See jayenne/regression/nightly_cmake_script.sh or
#   nightly_regression.csh
# - Summary: The script must do something like this:
#   [export work_dir=/full/path/to/working/dir]
#   ctest [-V] [-VV] -S /path/to/this/script.cmake,\
#     [Experimental|Nightly|Continuous],\
#     [Debug[,Coverage]|Release|RelWithDebInfo]

set( CTEST_PROJECT_NAME "Draco" )
message("source ${CTEST_SCRIPT_DIRECTORY}/draco_regression_macros.cmake" )
include( "${CTEST_SCRIPT_DIRECTORY}/draco_regression_macros.cmake" )
set_defaults()
parse_args()
find_tools()
set_svn_command("draco/trunk")

# Platform customization:
# 1. Ceilito - set TOOCHAIN_SETUP
platform_customization()

####################################################################
# The values in this section are optional you can either
# have them or leave them commented out
####################################################################

# this is the initial cache to use for the binary tree, be careful to escape
# any quotes inside of this string if you use it
set( CTEST_INITIAL_CACHE "
CMAKE_VERBOSE_MAKEFILE:BOOL=ON
CMAKE_BUILD_TYPE:STRING=${CTEST_BUILD_CONFIGURATION}
CMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
CTEST_CMAKE_GENERATOR:STRING=${CTEST_CMAKE_GENERATOR}
CTEST_USE_LAUNCHERS:STRING=${CTEST_USE_LAUNCHERS}
CTEST_TEST_TIMEOUT:STRING=${CTEST_TEST_TIMEOUT}

VENDOR_DIR:PATH=${VENDOR_DIR}
AUTODOCDIR:PATH=${AUTODOCDIR}
${TEST_PPE_BINDIR}
USE_CUDA:BOOL=${USE_CUDA}

${INIT_CACHE_PPE_PREFIX}
${TOOLCHAIN_SETUP}
${CT_CUSTOM_VARS}
# Set DRACO_DIAGNOSTICS and DRACO_TIMING:
${FULLDIAGNOSTICS}
${BOUNDS_CHECKING}
")

message("CTEST_INITIAL_CACHE =
----------------------------------------------------------------------
${CTEST_INITIAL_CACHE}
----------------------------------------------------------------------")

if( "${CTEST_CONFIGURE}x" STREQUAL "x" )
   set( CTEST_CONFIGURE OFF )
endif()
if( "${CTEST_BUILD}x" STREQUAL "x" )
   set( CTEST_BUILD OFF )
endif()
if( "${CTEST_TEST}x" STREQUAL "x" )
   set( CTEST_TEST OFF )
endif()
if( "${CTEST_SUBMIT}x" STREQUAL "x" )
   set( CTEST_SUBMIT OFF )
endif()
message("
--> Draco_Linux64.cmake modes:
    CTEST_CONFIGURE = ${CTEST_CONFIGURE}
    CTEST_BUILD     = ${CTEST_BUILD}
    CTEST_TEST      = ${CTEST_TEST}
    CTEST_SUBMIT    = ${CTEST_SUBMIT}
")

message("Parsing ${CTEST_SOURCE_DIRECTORY}/CTestCustom.cmake")
ctest_read_custom_files("${CTEST_SOURCE_DIRECTORY}")

if( "${CTEST_CONFIGURE}" STREQUAL "ON" )
  if( EXISTS "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt")
    # Empty the binary directory and recreate the CMakeCache.txt
    message( "ctest_empty_binary_directory( ${CTEST_BINARY_DIRECTORY} )" )
    ctest_empty_binary_directory( ${CTEST_BINARY_DIRECTORY} )
  endif()
  file( WRITE ${CTEST_BINARY_DIRECTORY}/CMakeCache.txt ${CTEST_INITIAL_CACHE} )
endif()

# Start
# if this is the 2nd (or 3rd) call to this script (i.e.: test and submit
# on different calls) then append the results.
# if( (${CTEST_SUBMIT} AND NOT ${CTEST_TEST}) OR
#     (${CTEST_TEST}   AND NOT ${CTEST_BUILD}) )
# if( ${CTEST_CONFIGURE} )
#   # Test and Submit on different calls -> The submit step should
#   # append the previous run.
#   message( "ctest_start( ${CTEST_MODEL} APPEND )")
#   ctest_start( ${CTEST_MODEL} APPEND )
# else()
#   message( "ctest_start( ${CTEST_MODEL} )")
#   ctest_start( ${CTEST_MODEL} )
# endif()

if( ${CTEST_CONFIGURE} )
  message( "ctest_start( ${CTEST_MODEL} )")
  ctest_start( ${CTEST_MODEL} )
else()
  message( "ctest_start( ${CTEST_MODEL} APPEND )")
  ctest_start( ${CTEST_MODEL} APPEND )
endif()

#if( NOT "${CTEST_MODEL}" MATCHES "Nightly" )
#   file  ( READ "${CTEST_BINARY_DIRECTORY}/Testing/TAG" tag_file )
#   string( REGEX MATCH "[^\n]*" BuildTag ${tag_file} )
#   set   ( TEST_OUTPUT_DIR "${CTEST_BINARY_DIRECTORY}/Testing/${BuildTag}" )
#   set   ( TEST_TEMP_DIR "${TEST_OUTPUT_DIR}/temp" )
#   file  ( MAKE_DIRECTORY ${TEST_TEMP_DIR} )
#endif()

# Update and Configure
if( "${CTEST_CONFIGURE}" STREQUAL "ON" )
   message( "ctest_update( SOURCE ${CTEST_SOURCE_DIRECTORY} RETURN_VALUE res )"  )
   ctest_update( SOURCE ${CTEST_SOURCE_DIRECTORY} RETURN_VALUE res )
   message( "Files updated: ${res}" )

   message( "setup_for_code_coverage()" )
   setup_for_code_coverage() # from draco_regression_macros.cmake
   message(  "ctest_configure()" )
   ctest_configure(
      BUILD        "${CTEST_BINARY_DIRECTORY}"
      RETURN_VALUE res) # LABELS label1 [label2]
endif()

# Build
if( "${CTEST_BUILD}" STREQUAL "ON" )
   # Autodoc
   if( "${CTEST_BUILD_CONFIGURATION}" STREQUAL "Release" AND
       "${CTEST_SITE}" MATCHES "ccscs7" AND
       NOT ${CMAKE_INSTALL_PREFIX} MATCHES "gcc-5.2.0" )
      message( "ctest_build( TARGET autodoc RETURN_VALUE res )" )
      ctest_build(
         TARGET autodoc
         RETURN_VALUE res
         NUMBER_ERRORS num_errors
         NUMBER_WARNINGS num_warnings
         )
      message( "build result:
   ${res}
   Build errors  : ${num_errors}
   Build warnings: ${num_warnings}" )
   endif()
   # Main build
   message( "ctest_build( TARGET install RETURN_VALUE res )" )
   ctest_build(
      TARGET install
      RETURN_VALUE res
      NUMBER_ERRORS num_errors
      NUMBER_WARNINGS num_warnings
      )
   message( "build result:
   ${res}
   Build errors  : ${num_errors}
   Build warnings: ${num_warnings}" )
endif()

# Test
if( "${CTEST_TEST}" STREQUAL "ON" )
#   if(ENABLE_C_CODECOVERAGE)
#     message( "ctest_test( SCHEDULE_RANDOM ON )" )
#     ctest_test( SCHEDULE_RANDOM ON )
#   else()
     find_num_procs_avail_for_running_tests() # returns num_test_procs
     message( "ctest_test( PARALLEL_LEVEL ${num_test_procs} SCHEDULE_RANDOM ON )" )
     ctest_test( PARALLEL_LEVEL ${num_test_procs}
                 SCHEDULE_RANDOM ON
                 # INCLUDE_LABEL <LABEL>
                 )
#   endif()

   # Process code coverage (bullseye) or dynamic analysis (valgrind)
   message("Processing code coverage or dynamic analysis")
   process_cc_or_da()
endif()

# Submit results
if( "${CTEST_SUBMIT}" STREQUAL "ON" )
   message( "ctest_submit()")
   ctest_submit()
endif()

message("end of ${CTEST_SCRIPT_NAME}.")
