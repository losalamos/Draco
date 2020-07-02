#-----------------------------*-cmake-*----------------------------------------#
# file   config/FindCaliper.cmake
# date   Tuesday, Jun 30, 2020, 08:11 am
# brief  Instructions for discovering the Caliper vendor libraries.
# note   Copyright (C) 2016-2020 Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

#.rst:
# FindCaliper
# ---------
#
# Find the LLNL Caliper library and header files.
#
# Caliper: A Performance Analysis Toolbox in a Library,
# https://github.com/llnl/Caliper
#
# Imported Targets
# ^^^^^^^^^^^^^^^^
#
# If Caliper is found, this module defines the following :prop_tgt:`IMPORTED`
# targets::
#
#  CALIPER::caliper        - The Caliper target
#
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# This module will set the following variables in your project::
#
#  CALIPER_FOUND          - True if Caliper found on the local system
#  CALIPER_INCLUDE_DIRS   - Location of Caliper header files.
#  CALIPER_LIBRARIES      - The Caliper libraries.
#  CALIPER_VERSION        - The version of the discovered Caliper install.
#
# Hints
# ^^^^^
#
# Set ``CALIPER_ROOT_DIR`` to a directory that contains a Caliper installation.
#
# This script expects to find libraries at ``$Caliper_ROOT_DIR/lib`` and the
# Caliper headers at ``$CALIPER_ROOT_DIR/include``.  The library directory may
# optionally provide Release and Debug folders.
#
# Cache Variables
# ^^^^^^^^^^^^^^^
#
# This module may set the following variables depending on platform and type of
# Caliper installation discovered.  These variables may optionally be set to
# help this module find the correct files::
#
#  CALIPER_LIBRARY        - Location of the Caliper library.
#  CALIPER_LIBRARY_DEBUG  - Location of the debug Caliper library (if any).
#
#------------------------------------------------------------------------------#

# Include these modules to handle the QUIETLY and REQUIRED arguments.
include(FindPackageHandleStandardArgs)

#=============================================================================
# If the user has provided ``CALIPER_ROOT_DIR``, use it!  Choose items found
# at this location over system locations.
if( EXISTS "$ENV{CALIPER_ROOT_DIR}" )
  file( TO_CMAKE_PATH "$ENV{CALIPER_ROOT_DIR}" CALIPER_ROOT_DIR )
  set( CALIPER_ROOT_DIR "${CALIPER_ROOT_DIR}" CACHE PATH
    "Prefix for Caliper installation." )
endif()

#=============================================================================
# Set CALIPER_INCLUDE_DIRS and CALIPER_LIBRARIES. Try to find the libraries at
# $CALIPER_ROOT_DIR (if provided) or in standard system locations.  These
# find_library and find_path calls will prefer custom locations over standard
# locations (HINTS).  If the requested file is not found at the HINTS location,
# standard system locations will be still be searched (/usr/lib64 (Redhat),
# lib/i386-linux-gnu (Debian)).

find_path( CALIPER_INCLUDE_DIR
  NAMES caliper/Caliper.h
  HINTS ${CALIPER_ROOT_DIR}/include
  PATH_SUFFIXES caliper
  )

set( CALIPER_LIBRARY_NAME caliper )

find_library(CALIPER_LIBRARY
  NAMES ${CALIPER_LIBRARY_NAME}
  PATHS ${CALIPER_ROOT_DIR}/lib64
  )
# Do we also have debug versions?
find_library( CALIPER_LIBRARY_DEBUG
  NAMES ${CALIPER_LIBRARY_NAME}
  HINTS ${CALIPER_ROOT_DIR}/lib64
  PATH_SUFFIXES Debug
)
set( CALIPER_INCLUDE_DIRS ${CALIPER_INCLUDE_DIR} )
set( CALIPER_LIBRARIES ${CALIPER_LIBRARY} )

#=============================================================================
# handle the QUIETLY and REQUIRED arguments and set CALIPER_FOUND to TRUE if
# all listed variables are TRUE.
find_package_handle_standard_args( Caliper
  FOUND_VAR
    CALIPER_FOUND
  REQUIRED_VARS
    CALIPER_INCLUDE_DIR
    CALIPER_LIBRARY
  VERSION_VAR
    CALIPER_VERSION
    )

mark_as_advanced( CALIPER_ROOT_DIR CALIPER_VERSION CALIPER_LIBRARY
  CALIPER_INCLUDE_DIR CALIPER_LIBRARY_DEBUG CALIPER_USE_PKGCONFIG
  CALIPER_CONFIG )

#=============================================================================
# Register imported libraries:
# 1. If we can find a Windows .dll file (or if we can find both Debug and
#    Release libraries), we will set appropriate target properties for these.
# 2. However, for most systems, we will only register the import location and
#    include directory.

# Look for dlls, or Release and Debug libraries.
if(WIN32)
  string( REPLACE ".lib" ".dll" CALIPER_LIBRARY_DLL
    "${CALIPER_LIBRARY}" )
  string( REPLACE ".lib" ".dll" CALIPER_LIBRARY_DEBUG_DLL
    "${CALIPER_LIBRARY_DEBUG}" )
endif()

if( CALIPER_FOUND AND NOT TARGET CALIPER::caliper )
  if( WIN32 )
    if( EXISTS "${CALIPER_LIBRARY_DLL}" )

      # Windows systems with dll libraries.
      add_library( CALIPER::caliper SHARED IMPORTED )

      # Windows with dlls, but only Release libraries.
      set_target_properties( CALIPER::caliper PROPERTIES
        IMPORTED_LOCATION_RELEASE         "${CALIPER_LIBRARY_DLL}"
        IMPORTED_IMPLIB                   "${CALIPER_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES     "${CALIPER_INCLUDE_DIRS}"
        IMPORTED_CONFIGURATIONS           Release
        IMPORTED_LINK_INTERFACE_LANGUAGES "C" )

      # If we have both Debug and Release libraries
      if( EXISTS "${CALIPER_LIBRARY_DEBUG_DLL}" )
        set_property( TARGET CALIPER::caliper APPEND PROPERTY
          IMPORTED_CONFIGURATIONS Debug )
        set_target_properties( CALIPER::caliper PROPERTIES
          IMPORTED_LOCATION_DEBUG           "${CALIPER_LIBRARY_DEBUG_DLL}"
          IMPORTED_IMPLIB_DEBUG             "${CALIPER_LIBRARY_DEBUG}" )
      endif()

    else()
      # Windows systems with static lib libraries.
      add_library( CALIPER::caliper SHARED IMPORTED )

      # Windows with dlls, but only Release libraries.
      set_target_properties( CALIPER::caliper PROPERTIES
        IMPORTED_LOCATION_RELEASE         "${CALIPER_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES     "${CALIPER_INCLUDE_DIRS}"
        IMPORTED_CONFIGURATIONS           Release
        IMPORTED_LINK_INTERFACE_LANGUAGES "C" )

      # If we have both Debug and Release libraries
      if( EXISTS "${CALIPER_LIBRARY_DEBUG}" )
        set_property( TARGET CALIPER::caliper APPEND PROPERTY
          IMPORTED_CONFIGURATIONS Debug )
        set_target_properties( CALIPER::caliper PROPERTIES
          IMPORTED_LOCATION_DEBUG           "${CALIPER_LIBRARY_DEBUG}" )
      endif()

    endif()

  else()

    # For all other environments (ones without dll libraries), create the
    # imported library targets.
    add_library( CALIPER::caliper    UNKNOWN IMPORTED )
    set_target_properties( CALIPER::caliper PROPERTIES
      IMPORTED_LOCATION                 "${CALIPER_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES     "${CALIPER_INCLUDE_DIRS}"
      IMPORTED_LINK_INTERFACE_LANGUAGES "C" )
  endif()
endif()

#------------------------------------------------------------------------------#
# End FindCaliper.cmake
#------------------------------------------------------------------------------#
