#-----------------------------*-cmake-*----------------------------------------#
# file   draco/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 April 28
# brief  Instructions for building root level Makefile.
# note   Copyright © 2010 Los Alamos National Security, All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# How to build Draco using CMake:
#------------------------------------------------------------------------------#
# Background:
# 
# There are 3 directories: source; build; target.
#   source - this is the location of cvs co files
#   build  - this is the location for makefiles, object files, etc.
#   target - this is the install location for final headers, libraries
#            and binaries. 
#
# Environment:
#
# - cmake must be in your path
# - The compiler is chosen by setting $CXX, $CC, $F90, etc.  On the
#   ICN this is done for you automatcially.  If not set, CMake will
#   use whatever it can find in your path (e.g.: /usr/bin/g++)
# - Vendor locations can be specified by setting $<VENDOR>_LIB_DIR and
#   $<VENDOR>_INC_DIR. For example, $MPI_INC_DIR, $GSL_LIB_DIR,
#   etc. When using the Draco modules feature, this is done
#   automatically. 
#
# Requirements:
#
# 1. Draco must find GSL and LAPACK.  If these are not found, the
#    configuration process will abort.
#
# Building:
# 
# 1. Change to your build directory.
# 2. [optional] Copy the $SRC_DIR/CMakeCache.txt file to the build
#    directory and edit the file as needed.
# 3. Run 'cmake $SRC_DIR'.  This generates Makefiles in the build
#    directory. 
# 4. [optional] Use command line arguments to cmake to configure your
#    build (instead of using CMakeCache.txt).  For example:
#    cmake ~/draco \
#          -DCMAKE_INSTALL_PREFIX=/var/tmp/draco-install-dir \
#          -DMPI_INC_DIR=/ccs/codes/mpi/openmpi/Linux64/1.3.3/include \
#          -DMPI_LIB_DIR=/ccs/codes/mpi/openmpi/Linux64/1.3.3/lib \
#          -DDRACO_DBC_LEVEL=0 \
#          -DCMAKE_BUILD_TYPE=Release
# 5. [optional] Run 'make edit_cache' to review your build setup.
#    Edit as needed. 
# 6. Run 'make -j16' to compile the code. (also try 'make VERBOSE=1')
# 7. Run 'ctest -j16' to run the tests.
# 8. Run 'make install' to install libraries in the $CMAKE_PREFIX_DIR.
# 9. [optional] Run 'make help' for a list of targets.
#------------------------------------------------------------------------------#

cmake_minimum_required( VERSION 2.8 FATAL_ERROR )
project( Draco CXX )

# Build system configuration files are located here.
set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config )

#
# The Draco version number.
#
set(DRACO_VERSION_MAJOR 6)
set(DRACO_VERSION_MINOR 0)
# set(DRACO_VERSION_PATCH 0) # only set this for a release
include(dracoVersion)

#
# Unit Test Setup
#
include(dracoTesting) 
# Top level test
add_test( SystemInformationNew 
	"${CMAKE_COMMAND}" --system-information  -G "${CMAKE_GENERATOR}" )

# Build options (static or shared, dbc-level
#   Insist() assertions    : always on
#   Require() preconditions: add +1 to DBC_LEVEL
#   Check() assertions     : add +2 to DBC_LEVEL
#   Ensure() postconditions: add +4 to DBC_LEVEL
set( DRACO_DBC_LEVEL "7" CACHE STRING "Design-by-Contract (0-7)?" )

#
# Packaging (CPack) options
# 
#TODO: Not a packaged product, so leave out for now:
#include( cpack_env )
# For example see CMake source dist: CMakeCPack.cmake.

#
# Build Draco components:
# 
add_subdirectory( src )
add_subdirectory( config )

#
# Dashboard/Regression setup
#
# CTest/CDash custom configuration file
configure_file(
  ${PROJECT_SOURCE_DIR}/CTestCustom.cmake
  ${PROJECT_BINARY_DIR}/CTestCustom.cmake
  COPYONLY )

#
# Extra custom commands
#
# alias for 'make configure'
add_custom_target( configure
   COMMAND ${CMAKE_MAKE_PROGRAM} rebuild_cache )   

# Doxygen documentation
option( BUILD_AUTODOC OFF )
if( BUILD_AUTODOC )
   add_subdirectory( autodoc )
endif()

# provide the cache file for reference
install( FILES ${PROJECT_BINARY_DIR}/CMakeCache.txt 
         DESTINATION config 
         RENAME Draco-CMakeCache.txt )
# install top level documents
install( 
   FILES 
      ChangeLog
      Copyright
      README.draco
      ToDo
   DESTINATION ${CMAKE_INSTALL_PREFIX}
)

# Add support for the 'uninstall' target
configure_file(
   "${CMAKE_CURRENT_SOURCE_DIR}/config/cmake_uninstall.cmake.in"
   "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
   @ONLY )
add_custom_target( uninstall
   "${CMAKE_COMMAND}" -P
   "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake" )

##---------------------------------------------------------------------------##
## Notes:

# 1. Use 'cmake -LAH' to list all user set variables.
# 2. Use 'make rebuild_cache' to regenerate config.h files, directories,
#   etc.
# 3. Use 'make edit_cache' to change user set variables.
# 4. Use 'make -jN' to build in parallel with N processes.
# 5. Use 'make test' or 'ctest -jN' to run the unit tests.
# 6. Use 'make install' to strip libraries and copy libraries and
#    headers to the INSTALL_PREFIX directory.


##---------------------------------------------------------------------------##
## End of CMakeLists.txt
##---------------------------------------------------------------------------##
