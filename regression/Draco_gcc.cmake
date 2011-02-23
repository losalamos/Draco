#
# Linux 64-bit
# G++/GCC/Gfortran 4.X.X
# Ref: http://www.cmake.org/Wiki/CMake_Scripting_Of_CTest

# 
# [export work_dir=/full/path/to/working/dir]
# ctest [-V] [-VV] -S /path/to/this/script.cmake,\
# [Experimental|Nightly|Continuous],[Debug[,Coverage]|Release|RelWithDebInfo]
#

set( CTEST_PROJECT_NAME "Draco" )
include( "${CTEST_SCRIPT_DIRECTORY}/draco_regression_macros.cmake" )

# Check inputs. Set default values for
#     CTEST_SOURCE_DIRECTORY
#     CTEST_BINARY_DIRECTORY
#     CMAKE_INSTALL_PREFIX
#     CMAKE_GENERATOR
#     CTEST_MODEL
#     CTEST_BUILD_CONFIGURATION
#     CTEST_START_WITH_EMPTY_BINARY_DIRECTORY
#     CTEST_CONTINUOUS_DURATION
#     CTEST_CONTINUOUS_MINIMUM_INTERVAL
#     VENDOR_DIR
#     CMAKE_GENERATOR
#     CTEST_NIGHTLY_START_TIME
#     CTEST_DROP_METHOD = http
#     CTEST_DROP_SITE   = coder.lanl.gov
#     CTEST_DROP_LOCATION = "/cdash/submit.php?project=${CTEST_PROJECT_NAME}" 
#     CTEST_DROP_SITE_CDASH = TRUE 
#     CTEST_CURL_OPTIONS    = CURLOPT_SSL_VERIFYPEER_OFF
#     sitename
set_defaults() # QUIET

# Based on command line, update values for
#     CTEST_MODEL
#     CTEST_BUILD_CONFIGURATION
#     build_name
#     enable_coverage
parse_args() # QUIET

# Finds tools and sets:
#     CTEST_CMD
#     CTEST_CVS_COMMAND
#     CTEST_CMAKE_COMMAND
find_tools() # QUIET

set( CTEST_CVS_CHECKOUT
  "${CTEST_CVS_COMMAND} -d /ccs/codes/radtran/cvsroot co -P -d source draco" )
#set( CTEST_CVS_CHECKOUT
#  "${CTEST_CVS_COMMAND} -d $ENV{USERNAME}@ccscs8.lanl.gov/ccs/codes/radtran/cvsroot co -P -d source draco" )
# under windows, consider: file:///z:/radiative/...

# set( CTEST_NOTES_FILE "path/to/file1" "/path/to/file2" )

####################################################################
# The values in this section are optional you can either
# have them or leave them commented out
####################################################################

# Test Coverage setup
if( ENABLE_C_CODECOVERAGE )
   set( ENV{COVFILE} "${CTEST_BINARY_DIRECTORY}/CMake.cov" )
   set( ENV{COVDIRCFG} "${CTEST_SCRIPT_DIRECTORY}/covclass_cmake.cfg" )
   set( ENV{COVFNCFG} "${CTEST_SCRIPT_DIRECTORY}/covclass_cmake.cfg" )
   set( ENV{COVCLASSCFG} "${CTEST_SCRIPT_DIRECTORY}/covclass_cmake.cfg" )
   set( ENV{COVSRCCFG} "${CTEST_SCRIPT_DIRECTORY}/covclass_cmake.cfg" )
endif()

# this is the initial cache to use for the binary tree, be careful to escape
# any quotes inside of this string if you use it
set( CTEST_INITIAL_CACHE "
CMAKE_VERBOSE_MAKEFILE:BOOL=ON
CMAKE_BUILD_TYPE:STRING=${CTEST_BUILD_CONFIGURATION}
CMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
CTEST_CMAKE_GENERATOR:STRING=${CTEST_CMAKE_GENERATOR}
CTEST_USE_LAUNCHERS:STRING=${CTEST_USE_LAUNCHERS}

ENABLE_C_CODECOVERAGE:BOOL=${ENABLE_C_CODECOVERAGE}
ENABLE_Fortran_CODECOVERAGE:BOOL=${ENABLE_Fortran_CODECOVERAGE}
VENDOR_DIR:PATH=/ccs/codes/radtran/vendors/Linux64
")
# set( CTEST_INITIAL_CACHE "
# BUILDNAME:STRING=${build_name}
# CMAKE_MAKE_PROGRAM:FILEPATH=${MAKECOMMAND}
# CVSCOMMAND:FILEPATH=${CTEST_CVS_COMMAND}
# ENABLE_C_CODECOVERAGE:BOOL=${ENABLE_C_CODECOVERAGE}
# ENABLE_Fortran_CODECOVERAGE:BOOL=${ENABLE_Fortran_CODECOVERAGE}
# MAKECOMMAND:FILEPATH=${MAKECOMMAND}
# SITE:STRING=${sitename}
# SVNCOMMAND:FILEPATH=${CTEST_CVS_COMMAND}
# VENDOR_DIR:PATH=/ccs/codes/radtran/vendors/Linux64
# MEMORYCHECK_COMMAND:FILEPATH=${MEMORYCHECK_COMMAND}
# MEMORYCHECK_COMMAND_OPTIONS:STRING=${MEMORYCHECK_COMMAND_OPTIONS}
# MEMORYCHECK_SUPPRESSIONS_FILE:FILEPATH=${MEMORYCHECK_SUPPRESSIONS_FILE}
# CC:FILEPATH=${CC}
# CXX:FILEPATH=${CXX}

# ")
message("CTEST_INITIAL_CACHE =  
${CTEST_INITIAL_CACHE}
")

ctest_empty_binary_directory( ${CTEST_BINARY_DIRECTORY} )
file( WRITE ${CTEST_BINARY_DIRECTORY}/CMakeCache.txt ${CTEST_INITIAL_CACHE} )

# set any extra environment variables to use during the execution of
# the script here: 
#set( ENV{FC} $ENV{F90} )
set( VERBOSE ON )
set( CTEST_OUTPUT_ON_FAILURE ON )

# Set the CTEST_COMMAND
setup_ctest_commands() # QUIET

# Install the files
message( STATUS "Installing files to ${CMAKE_INSTALL_PREFIX}..." )
execute_process( 
   COMMAND           ${CMAKE_MAKE_PROGRAM} install
   WORKING_DIRECTORY ${CTEST_BINARY_DIRECTORY}
   )

message("end of ${CTEST_SCRIPT_NAME}.")

