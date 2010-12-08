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
  "${CTEST_CVS_COMMAND} -d $ENV{USERNAME}@ccscs8.lanl.gov/ccs/codes/radtran/cvsroot co -P -d source draco" )
# under windows, consider: file:///z:/radiative/...

# Set the CTEST_COMMAND
setup_ctest_commands() # QUIET

####################################################################
# The values in this section are optional you can either
# have them or leave them commented out
####################################################################

# this is the initial cache to use for the binary tree, be careful to escape
# any quotes inside of this string if you use it
set( CTEST_INITIAL_CACHE "
DRACO_BUILD_TESTS:BOOL=ON
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

# Ideas to try

#SET (CTEST_INITIAL_CACHE "
#CMAKE_GENERATOR:INTERNAL=Visual Studio 8 2005
#CMAKE_MAKE_PROGRAM:FILEPATH=C:/Program Files/Microsoft Visual Studio 8/Common7/IDE/devenv.com
#MAKECOMMAND:STRING=\"C:/Program Files/Microsoft Visual Studio 8/Common7/IDE/devenv.com\" VTK.sln /build Release /project ALL_BUILD

