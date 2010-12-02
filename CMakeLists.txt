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

#
# Build the packages:
# 
add_subdirectory( src )
add_subdirectory( config )

# Top level tests
add_test( SystemInformationNew 
	"${CMAKE_COMMAND}" --system-information  -G "${CMAKE_GENERATOR}" )

# CTest/CDash custom configuration file
#configure_file(
#  ${PROJECT_SOURCE_DIR}/CTestCustom.cmake
#  ${PROJECT_BINARY_DIR}/CTestCustom.cmake
#  COPYONLY
#  )

#TODO: Leave out for now:
#
# Provide some documentation
# 
# install( 
  # FILES
    # ${PROJECT_SOURCE_DIR}/ReleaseNote.txt 
    # ${PROJECT_SOURCE_DIR}/Copyright.txt 
  # DESTINATION .
  # COMPONENT docs
  #OPTIONAL
  # )

# Define release targets (Win32 only)
#include( release )

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

# set PATH=C:\Program Files (x86)\Microsoft Visual Studio 10.0\Common7\IDE\;C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\BIN;C:\Program Files (x86)\Microsoft Visual Studio 10.0\Common7\Tools;C:\Windows\Microsoft.NET\Framework\v4.0.30319;C:\Windows\Microsoft.NET\Framework\v3.5;C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\VCPackages;C:\Program Files (x86)\Microsoft SDKs\Windows\v7.0A\bin\NETFX 4.0 Tools;C:\Program Files (x86)\Microsoft SDKs\Windows\v7.0A\bin;C:\Program Files\Common Files\Microsoft Shared\Windows Live;C:\Program Files (x86)\Common Files\Intel\Shared Files\fortran\bin\Intel64;C:\Windows\system32;C:\Windows;C:\Windows\System32\Wbem;C:\Program Files (x86)\QT Lite\QTSystem;c:\Program Files (x86)\Microsoft SQL Server\90\Tools\binn\;c:\Program Files (x86)\Microsoft SQL Server\100\Tools\Binn\;c:\Program Files\Microsoft SQL Server\100\Tools\Binn\;c:\Program Files\Microsoft SQL Server\100\DTS\Binn;C:\Program Files (x86)\CMake 2.8\bin;C:\Program Files (x86)\QT Lite\QTSystem\;C:\Program Files (x86)\MPICH2\bin

# Notes:

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
