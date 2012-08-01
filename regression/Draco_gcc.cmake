#
# Linux 64-bit
# Ref: http://www.cmake.org/Wiki/CMake_Scripting_Of_CTest

cmake_minimum_required(VERSION 2.6.3 FATAL_ERROR)

# Use:
# - See jayenne/regression/nightly_cmake_script.sh or
#   nightly_regression.csh
# - Summary: The script must do something like this:
#   [export work_dir=/full/path/to/working/dir]
#   ctest [-V] [-VV] -S /path/to/this/script.cmake,\
#     [Experimental|Nightly|Continuous],\
#     [Debug[,Coverage]|Release|RelWithDebInfo]

set( CTEST_PROJECT_NAME "Draco" )
include( "${CTEST_SCRIPT_DIRECTORY}/draco_regression_macros.cmake" )
set_defaults()
parse_args() 
find_tools()
set_svn_command("draco/trunk")

# Platform customization:
# 1. Ceilito - set TOOCHAIN_SETUP
# 2. Roadrunner - set TOOCHAIN_SETUP and INIT_CACHE_PPE_PREFIX
platform_customization()
if( "${sitename}" MATCHES "RoadRunner" )
   if( NOT "$ENV{CXX}" MATCHES "ppu-g[+][+]" )
      string( REPLACE "${CTEST_MODEL}_g++" "${CTEST_MODEL}_ppu-g++" TEST_PPE_BINDIR $ENV{work_dir} )
      set(TEST_PPE_BINDIR "TEST_PPE_BINDIR:PATH=${TEST_PPE_BINDIR}/build/src/device/test")
   endif()
endif()

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

ENABLE_C_CODECOVERAGE:BOOL=${ENABLE_C_CODECOVERAGE}
ENABLE_Fortran_CODECOVERAGE:BOOL=${ENABLE_Fortran_CODECOVERAGE}
VENDOR_DIR:PATH=/ccs/codes/radtran/vendors/Linux64
${TEST_PPE_BINDIR}

${INIT_CACHE_PPE_PREFIX}
${TOOLCHAIN_SETUP}
${CT_CUSTOM_VARS}
")

message("CTEST_INITIAL_CACHE =  
----------------------------------------------------------------------
${CTEST_INITIAL_CACHE}
----------------------------------------------------------------------")

message("Parsing ${CTEST_SOURCE_DIRECTORY}/CTestCustom.cmake")
ctest_read_custom_files("${CTEST_SOURCE_DIRECTORY}")

if( "${CTEST_SUBMITONLY}" STREQUAL "ON" )
   # Start
   message( "ctest_start( ${CTEST_MODEL} )")
   ctest_start( ${CTEST_MODEL} )

   # Submit results
   message( "ctest_submit()")
   ctest_submit()

else( "${CTEST_SUBMITONLY}" STREQUAL "ON" )

   # Empty the binary directory and recreate the CMakeCache.txt
   message( "ctest_empty_binary_directory( ${CTEST_BINARY_DIRECTORY} )" )
   ctest_empty_binary_directory( ${CTEST_BINARY_DIRECTORY} )
   file( WRITE ${CTEST_BINARY_DIRECTORY}/CMakeCache.txt ${CTEST_INITIAL_CACHE} )

   # Start
   message( "ctest_start( ${CTEST_MODEL} )")
   ctest_start( ${CTEST_MODEL} )

   # Update
   message( "ctest_update( SOURCE ${CTEST_SOURCE_DIRECTORY} RETURN_VALUE res )"  )
   ctest_update( SOURCE ${CTEST_SOURCE_DIRECTORY} RETURN_VALUE res )
   message( "Files updated: ${res}" )

   # Configure
   message( "setup_for_code_coverage()" )
   setup_for_code_coverage() # from draco_regression_macros.cmake
   message(  "ctest_configure()" )
   ctest_configure(
      BUILD        "${CTEST_BINARY_DIRECTORY}"
      RETURN_VALUE res) # LABELS label1 [label2] 

   # Build
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

   # Test
   message( "ctest_test( PARALLEL_LEVEL ${MPIEXEC_MAX_NUMPROCS} SCHEDULE_RANDOM ON )" )
   ctest_test( PARALLEL_LEVEL ${MPIEXEC_MAX_NUMPROCS} SCHEDULE_RANDOM ON ) 

   # Process code coverage (bullseye) or dynamic analysis (valgrind)
   message("Processing code coverage or dynamic analysis")
   process_cc_or_da()

   if( NOT "${CTEST_NOSUBMIT}" STREQUAL "ON" )
      # Submit results
      message( "ctest_submit()")
      ctest_submit()
   endif()

   # Install the files
   # message(  "Installing files to ${CMAKE_INSTALL_PREFIX}..." )
   # execute_process( 
   #    COMMAND           ${CMAKE_MAKE_PROGRAM} install
   #    WORKING_DIRECTORY ${CTEST_BINARY_DIRECTORY}
   #    )
endif( "${CTEST_SUBMITONLY}" STREQUAL "ON" )

message("end of ${CTEST_SCRIPT_NAME}.")

