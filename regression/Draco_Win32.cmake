#-----------------------------*-cmake-*----------------------------------------#
# file   draco/regression/Draco_Win32.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   2016 June 14
# brief  CTest regression script for Draco on Win32.
# note   Copyright (C) 2016-2019 Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# Ref: http://www.cmake.org/Wiki/CMake_Scripting_Of_CTest

cmake_minimum_required(VERSION 3.9.0)

# Use:
# - See draco/regression/win32-regression-master.bat
# - Summary: The script must do something like this:
#   [export work_dir=/full/path/to/working/dir]
#   ctest [-V] [-VV] -S /path/to/this/script.cmake,\
#     [Experimental|Nightly|Continuous],\
#     [Debug[,Coverage|,DynamicAnalysis]|Release|RelWithDebInfo]

set( CTEST_PROJECT_NAME "Draco" )
message("source ${CTEST_SCRIPT_DIRECTORY}/draco_regression_macros.cmake" )
include( "${CTEST_SCRIPT_DIRECTORY}/draco_regression_macros.cmake" )
set_defaults()
parse_args()
find_tools()
set_git_command("Draco.git")

# Make machine name lower case
string( TOLOWER "${CTEST_SITE}" CTEST_SITE )

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
# CMAKE_MAKE_PROGRAM:FILEPATH=${MAKECOMMAND}
${TEST_PPE_BINDIR}
WITH_CUDA:BOOL=${WITH_CUDA}

${INIT_CACHE_PPE_PREFIX}
${TOOLCHAIN_SETUP}
# Set DRACO_DIAGNOSTICS and DRACO_TIMING:
${FULLDIAGNOSTICS}
${BOUNDS_CHECKING}
")

message("CTEST_INITIAL_CACHE =
----------------------------------------------------------------------
${CTEST_INITIAL_CACHE}
----------------------------------------------------------------------")

message("
--> Draco_Linux64.cmake modes:
    CTEST_CONFIGURE = ${CTEST_CONFIGURE}
    CTEST_BUILD     = ${CTEST_BUILD}
    CTEST_TEST      = ${CTEST_TEST}
    CTEST_SUBMIT    = ${CTEST_SUBMIT}
")

message("Parsing ${CTEST_SOURCE_DIRECTORY}/CTestCustom.cmake")
ctest_read_custom_files("${CTEST_SOURCE_DIRECTORY}")

# Init
if( ${CTEST_CONFIGURE} )
  if( EXISTS "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt")
    # Empty the binary directory and recreate the CMakeCache.txt
    message( "ctest_empty_binary_directory( ${CTEST_BINARY_DIRECTORY} )" )
    ctest_empty_binary_directory( ${CTEST_BINARY_DIRECTORY} )
  endif()
  # dummy command to give the file system time to catch up before creating
  # CMakeCache.txt.
  # file( WRITE $ENV{TEMP}/foo.txt ${CTEST_INITIAL_CACHE} )
  file( WRITE ${CTEST_BINARY_DIRECTORY}/CMakeCache.txt ${CTEST_INITIAL_CACHE} )
endif()

if( ${CTEST_CONFIGURE} )
  message( "ctest_start( ${CTEST_MODEL} )")
  ctest_start( ${CTEST_MODEL} )
else()
  message( "ctest_start( ${CTEST_MODEL} APPEND )")
  ctest_start( ${CTEST_MODEL} APPEND )
endif()

# Update and Configure
if( ${CTEST_CONFIGURE} )
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

# Autodoc
if( ${CTEST_BUILD} AND ${CTEST_AUTODOC} )
  message( "ctest_build(
   TARGET autodoc
   NUMBER_ERRORS num_errors
   NUMBER_WARNINGS num_warnings
   FLAGS -j 24
   RETURN_VALUE res )" )
  ctest_build(
    TARGET autodoc
    RETURN_VALUE res
    NUMBER_ERRORS num_errors
    NUMBER_WARNINGS num_warnings
    FLAGS "-j 24"
    )
  message( "build result:
   ${res}
   Build errors  : ${num_errors}
   Build warnings: ${num_warnings}" )
endif()

# Build
if( ${CTEST_BUILD} )
   message( "ctest_build(
   TARGET install
   RETURN_VALUE res
   NUMBER_ERRORS num_errors
   NUMBER_WARNINGS num_warnings )" )
   ctest_build(
      TARGET install
      RETURN_VALUE res
      NUMBER_ERRORS num_errors
      NUMBER_WARNINGS num_warnings )
   message( "build result:
   ${res}
   Build errors  : ${num_errors}
   Build warnings: ${num_warnings}" )
endif()

# Test
if( ${CTEST_TEST} )

  find_num_procs_avail_for_running_tests() # returns num_test_procs
  set( ctest_test_options "SCHEDULE_RANDOM ON" )
  string( APPEND ctest_test_options " PARALLEL_LEVEL ${num_test_procs}" )

  # if we are running on a machine that openly shares resources, use the
  # TEST_LOAD feature to limit the number of cores used while testing. For
  # machines that run schedulers, the whole allocation is available so there is
  # no need to limit the load.
  if( "${CTEST_SITE}" MATCHES "ccscs" )
    string( APPEND ctest_test_options " TEST_LOAD ${max_system_load}" )
  endif()

  message( "ctest_test( ${ctest_test_options} )" )
  # convert string to a cmake list
  string( REPLACE " " ";" ctest_test_options "${ctest_test_options}" )
  ctest_test( ${ctest_test_options} )

  # Process code coverage (bullseye) or dynamic analysis (valgrind)
  message("Processing code coverage or dynamic analysis")
  process_cc_or_da()
endif()

# Submit results
if( ${CTEST_SUBMIT} )
   message( "ctest_submit()")
   ctest_submit()
endif()

message("end of ${CTEST_SCRIPT_NAME}.")

#------------------------------------------------------------------------------#
# End Draco_Win32.cmake
#------------------------------------------------------------------------------#
