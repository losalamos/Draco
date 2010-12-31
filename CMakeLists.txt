#-----------------------------*-cmake-*----------------------------------------#
# file   draco/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 April 28
# brief  Instructions for building root level Makefile.
# note   Copyright © 2010 Los Alamos National Security, All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

cmake_minimum_required( VERSION 2.8 FATAL_ERROR )
project( Draco CXX )

# Build system configuration files are located here.
set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config )

# The Draco version number.
include(dracoVersion)

# Select packages to build.
option( DRACO_BUILD_JAYENNE_COMPONENTS ON )
option( DRACO_BUILD_CAPSAICIN_COMPONENTS ON )

# add uninstall target
#include( uninstall )
#provide_uninstall_target()

# Unit Test Setup
option( DRACO_BUILD_TESTS "Should we compile the tests?" ON )
if( DRACO_BUILD_TESTS )
  # include( dracoTesting.cmake )
  include(CTest)
  enable_testing()
  if( "${DRACO_NUM_CORES}none" STREQUAL "none"  )
     add_custom_target( check
        COMMAND ${CMAKE_MAKE_COMMAND} test )   
  else()
     add_custom_target( check
        COMMAND ${CMAKE_CTEST_COMMAND} -j ${DRACO_NUM_CORES} )   
  endif()
endif()   

# Build options (static or shared, dbc-level
# Insist() assertions    : always on
# Require() preconditions: add +1 to DBC_LEVEL
# Check() assertions     : add +2 to DBC_LEVEL
# Ensure() postconditions: add +4 to DBC_LEVEL
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

# Top level tests
#add_test( SystemInformationNew 
#	"${CMAKE_COMMAND}" --system-information  -G "${CMAKE_GENERATOR}" )

# CTest/CDash custom configuration file
configure_file(
  ${PROJECT_SOURCE_DIR}/CTestCustom.cmake
  ${PROJECT_BINARY_DIR}/CTestCustom.cmake
  COPYONLY
  )

# Extra custom commands

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
