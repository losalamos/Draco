#
# Linux 64-bit
# G++/GCC/Gfortran 4.X.X
# Ref: http://www.cmake.org/Wiki/CMake_Scripting_Of_CTest

# 
# [export work_dir=/full/path/to/working/dir]
# ctest [-V] [-VV] -S /path/to/this/script.cmake,\
# [Experimental|Nightly|Continuous],[Debug[,Coverage]|Release|RelWithDebInfo]
#

include( "${CTEST_SCRIPT_DIRECTORY}/draco_regression_macros.cmake" )

# Check inputs. Set default values for
#     CTEST_SOURCE_DIRECTORY
#     CTEST_BINARY_DIRECTORY
#     CMAKE_INSTALL_PREFIX
#     CMAKE_GENERATOR
#     dashboard_type
#     build_type
#     CTEST_START_WITH_EMPTY_BINARY_DIRECTORY
#     CTEST_CONTINUOUS_DURATION
#     CTEST_CONTINUOUS_MINIMUM_INTERVAL
#     VENDOR_DIR
#     CMAKE_GENERATOR
#     sitename
set_defaults() # QUIET

# Based on command line, update values for
#     dashboard_type
#     build_type
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

# Set the CTEST_COMMAND
setup_ctest_commands() # QUIET

####################################################################
# The values in this section are optional you can either
# have them or leave them commented out
####################################################################

# clear the binary directory and create an initial cache
ctest_empty_binary_directory( ${CTEST_BINARY_DIRECTORY} )

# Test Coverage setup
if( ENABLE_C_CODECOVERAGE )
   find_program( COV01 cov01 )
   get_filename_component( beyedir ${COV01} PATH )
   set( CC ${beyedir}/gcc )
   set( CXX ${beyedir}/g++ )
   set( ENV{CC} ${beyedir}/gcc )
   set( ENV{CXX} ${beyedir}/g++ )

   # Set the coverage data file.
   set( ENV{COVFILE} "${CTEST_BINARY_DIRECTORY}/CMake.cov" )
   set( ENV{COVDIRCFG} "${CTEST_SCRIPT_DIRECTORY}/covclass_cmake.cfg" )
   set( ENV{COVFNCFG} "${CTEST_SCRIPT_DIRECTORY}/covclass_cmake.cfg" )
   set( ENV{COVCLASSCFG} "${CTEST_SCRIPT_DIRECTORY}/covclass_cmake.cfg" )
   set( ENV{COVSRCCFG} "${CTEST_SCRIPT_DIRECTORY}/covclass_cmake.cfg" )

   # turn off coverage for configure step
   set( RES 1 )
   execute_process(COMMAND ${COV01} -1 RESULT_VARIABLE RES)
   if( RES )
      message(FATAL_ERROR "could not run cov01 -1")
   else()
      message(STATUS "BullseyeCoverage turned on")
   endif()
  
endif()

# this is the initial cache to use for the binary tree, be careful to escape
# any quotes inside of this string if you use it
set( CTEST_INITIAL_CACHE "
BUILD_TESTING:BOOL=ON
VERBOSE:BOOL=ON

BUILDNAME:STRING=${build_name}
CMAKE_BUILD_TYPE:STRING=${build_type}
CMAKE_GENERATOR:STRING=${CMAKE_GENERATOR}
CMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
CMAKE_MAKE_PROGRAM:FILEPATH=${MAKECOMMAND}
CVSCOMMAND:FILEPATH=${CTEST_CVS_COMMAND}
ENABLE_C_CODECOVERAGE:BOOL=${ENABLE_C_CODECOVERAGE}
ENABLE_Fortran_CODECOVERAGE:BOOL=${ENABLE_Fortran_CODECOVERAGE}
MAKECOMMAND:FILEPATH=${MAKECOMMAND} -j8
SITE:STRING=${sitename}
SVNCOMMAND:FILEPATH=${CTEST_CVS_COMMAND}
VENDOR_DIR:PATH=/ccs/codes/radtran/vendors/Linux64
MEMORYCHECK_COMMAND:FILEPATH=${MEMORYCHECK_COMMAND}
CC:FILEPATH=${CC}
CXX:FILEPATH=${CXX}
")
message("

CTEST_INITIAL_CACHE = ${CTEST_INITIAL_CACHE}
")

# set any extra environment variables to use during the execution of
# the script here: 
set( CTEST_ENVIRONMENT
  FC=$ENV{F90}
  VERBOSE=ON
  CTEST_OUTPUT_ON_FAILURE=ON
)

message("end of ${CTEST_SCRIPT_NAME}.")


