#-----------------------------*-cmake-*----------------------------------------#
# file   config/FindEOSPAC.cmake
# date   2017 February 28
# brief  Instructions for discovering the EOSPAC vendor libraries.
# note   Copyright (C) 2016-2019 Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

#.rst:
# FindEOSPAC
# ---------
#
# Find the LANL EOSPAC library and header files.
#
# A collection of C routines that can be used to access the Sesame data
# library. https://laws.lanl.gov/projects/data/eos.html
#
# Imported Targets
# ^^^^^^^^^^^^^^^^
#
# If EOSPAC is found, this module defines the following :prop_tgt:`IMPORTED`
# targets::
#
#  EOSPAC::eospac        - The EOSPAC library.
#
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# This module will set the following variables in your project::
#
#  EOSPAC_FOUND          - True if EOSPAC found on the local system
#  EOSPAC_INCLUDE_DIRS   - Location of EOSPAC header files.
#  EOSPAC_LIBRARIES      - The EOSPAC libraries.
#  EOSPAC_VERSION        - The version of the discovered EOSPAC install.
#
# Hints
# ^^^^^
#
# Set ``EOSPAC_ROOT_DIR`` to a directory that contains a EOSPAC installation.
#
# This script expects to find libraries at ``$EOSPAC_ROOT_DIR/lib`` and the
# EOSPAC headers at ``$EOSPAC_ROOT_DIR/include``.  The library directory may
# optionally provide Release and Debug folders.
#
# Cache Variables
# ^^^^^^^^^^^^^^^
#
# This module may set the following variables depending on platform and type of
# EOSPAC installation discovered.  These variables may optionally be set to
# help this module find the correct files::
#
#  EOSPAC_LIBRARY        - Location of the EOSPAC library.
#  EOSPAC_LIBRARY_DEBUG  - Location of the debug EOSPAC library (if any).
#
#------------------------------------------------------------------------------#

# Include these modules to handle the QUIETLY and REQUIRED arguments.
include(FindPackageHandleStandardArgs)

#=============================================================================
# If the user has provided ``EOSPAC_ROOT_DIR``, use it!  Choose items found
# at this location over system locations.
if( EXISTS "$ENV{EOSPAC_ROOT_DIR}" )
  file( TO_CMAKE_PATH "$ENV{EOSPAC_ROOT_DIR}" EOSPAC_ROOT_DIR )
  set( EOSPAC_ROOT_DIR "${EOSPAC_ROOT_DIR}" CACHE PATH
    "Prefix for EOSPAC installation." )
endif()

#=============================================================================
# Set EOSPAC_INCLUDE_DIRS and EOSPAC_LIBRARIES. Try to find the libraries at
# $EOSPAC_ROOT_DIR (if provided) or in standard system locations.  These
# find_library and find_path calls will prefer custom locations over standard
# locations (HINTS).  If the requested file is not found at the HINTS location,
# standard system locations will be still be searched (/usr/lib64 (Redhat),
# lib/i386-linux-gnu (Debian)).

find_path( EOSPAC_INCLUDE_DIR
  NAMES eos_Interface.h
  HINTS ${EOSPAC_ROOT_DIR}/include
  PATH_SUFFIXES Release Debug
  )

# if (APPLE)
#     set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;.dylib")
# endif()

set( EOSPAC_LIBRARY_NAME eospac6)

find_library(EOSPAC_LIBRARY
  NAMES ${EOSPAC_LIBRARY_NAME}
  PATHS ${EOSPAC_ROOT_DIR}/lib
  PATH_SUFFIXES Release Debug
  )

# Do we also have debug versions?
find_library( EOSPAC_LIBRARY_DEBUG
  NAMES ${EOSPAC_LIBRARY_NAME}
  HINTS ${EOSPAC_ROOT_DIR}/lib 
  PATH_SUFFIXES Debug
)
set( EOSPAC_INCLUDE_DIRS ${EOSPAC_INCLUDE_DIR} )
set( EOSPAC_LIBRARIES ${EOSPAC_LIBRARY} )

# Try to find the version.
if( NOT EOSPAC_VERSION )
#   if( EXISTS "${EOSPAC_INCLUDE_DIRS}/eospac.h" )
#     file( STRINGS "${EOSPAC_INCLUDE_DIRS}/eospac.h" eospac_h_major
#         REGEX "define EOSPAC_VER_MAJOR" )
#     file( STRINGS "${EOSPAC_INCLUDE_DIRS}/eospac.h" eospac_h_minor
#         REGEX "define EOSPAC_VER_MINOR" )
#     file( STRINGS "${EOSPAC_INCLUDE_DIRS}/eospac.h" eospac_h_subminor
#         REGEX "define EOSPAC_VER_SUBMINOR" )
#     string( REGEX REPLACE ".*([0-9]+)" "\\1" EOSPAC_MAJOR ${eospac_h_major} )
#     string( REGEX REPLACE ".*([0-9]+)" "\\1" EOSPAC_MINOR ${eospac_h_minor} )
#     string( REGEX REPLACE ".*([0-9]+)" "\\1" EOSPAC_SUBMINOR ${eospac_h_subminor} )
#   endif()
  # We might also try scraping the directory name for a regex match
  # "eospac-X.X.X"
  if( NOT EOSPAC_MAJOR )
    string( REGEX REPLACE ".*eospac-([0-9]+).([0-9]+).([0-9]+).*" "\\1"
      EOSPAC_MAJOR ${EOSPAC_INCLUDE_DIR} )
    string( REGEX REPLACE ".*eospac-([0-9]+).([0-9]+).([0-9]+).*" "\\2"
      EOSPAC_MINOR ${EOSPAC_INCLUDE_DIR} )
    string( REGEX REPLACE ".*eospac-([0-9]+).([0-9]+).([0-9]+).*" "\\3"
      EOSPAC_SUBMINOR ${EOSPAC_INCLUDE_DIR} )
    set( EOSPAC_VERSION "${EOSPAC_MAJOR}.${EOSPAC_MINOR}.${EOSPAC_SUBMINOR}")
  endif()
  # Another option is `strings libeospac6.a | grep eos_version_name`
endif()

#=============================================================================
# handle the QUIETLY and REQUIRED arguments and set EOSPAC_FOUND to TRUE if
# all listed variables are TRUE.
find_package_handle_standard_args( EOSPAC
  FOUND_VAR
    EOSPAC_FOUND
  REQUIRED_VARS
    EOSPAC_INCLUDE_DIR
    EOSPAC_LIBRARY
  VERSION_VAR
    EOSPAC_VERSION
    )

mark_as_advanced( EOSPAC_ROOT_DIR EOSPAC_VERSION EOSPAC_LIBRARY
  EOSPAC_INCLUDE_DIR EOSPAC_LIBRARY_DEBUG EOSPAC_USE_PKGCONFIG EOSPAC_CONFIG )

#=============================================================================
# Register imported libraries:
# 1. If we can find a Windows .dll file (or if we can find both Debug and
#    Release libraries), we will set appropriate target properties for these.
# 2. However, for most systems, we will only register the import location and
#    include directory.

# Look for dlls, or Release and Debug libraries.
if(WIN32)
  string( REPLACE ".lib" ".dll" EOSPAC_LIBRARY_DLL
    "${EOSPAC_LIBRARY}" )
  string( REPLACE ".lib" ".dll" EOSPAC_LIBRARY_DEBUG_DLL
    "${EOSPAC_LIBRARY_DEBUG}" )
endif()

if( EOSPAC_FOUND AND NOT TARGET EOSPAC::eospac )
  if( WIN32 )
    if( EXISTS "${EOSPAC_LIBRARY_DLL}" )

      # Windows systems with dll libraries.
      add_library( EOSPAC::eospac SHARED IMPORTED )

      # Windows with dlls, but only Release libraries.
      set_target_properties( EOSPAC::eospac PROPERTIES
        IMPORTED_LOCATION_RELEASE         "${EOSPAC_LIBRARY_DLL}"
        IMPORTED_IMPLIB                   "${EOSPAC_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES     "${EOSPAC_INCLUDE_DIRS}"
        IMPORTED_CONFIGURATIONS           Release
        IMPORTED_LINK_INTERFACE_LANGUAGES "C" )

      # If we have both Debug and Release libraries
      if( EXISTS "${EOSPAC_LIBRARY_DEBUG_DLL}" )
        set_property( TARGET EOSPAC::eospac APPEND PROPERTY
          IMPORTED_CONFIGURATIONS Debug )
        set_target_properties( EOSPAC::eospac PROPERTIES
          IMPORTED_LOCATION_DEBUG           "${EOSPAC_LIBRARY_DEBUG_DLL}"
          IMPORTED_IMPLIB_DEBUG             "${EOSPAC_LIBRARY_DEBUG}" )
      endif()

    else()
      # Windows systems with static lib libraries.
      add_library( EOSPAC::eospac STATIC IMPORTED )

      # Windows with dlls, but only Release libraries.
      set_target_properties( EOSPAC::eospac PROPERTIES
        IMPORTED_LOCATION_RELEASE         "${EOSPAC_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES     "${EOSPAC_INCLUDE_DIRS}"
        IMPORTED_CONFIGURATIONS           Release
        IMPORTED_LINK_INTERFACE_LANGUAGES "C" )

      # If we have both Debug and Release libraries
      if( EXISTS "${EOSPAC_LIBRARY_DEBUG}" )
        set_property( TARGET EOSPAC::eospac APPEND PROPERTY
          IMPORTED_CONFIGURATIONS Debug )
        set_target_properties( EOSPAC::eospac PROPERTIES
          IMPORTED_LOCATION_DEBUG           "${EOSPAC_LIBRARY_DEBUG}" )
      endif()

    endif()
    
  else()

    # For all other environments (ones without dll libraries), create the
    # imported library targets.
    add_library( EOSPAC::eospac    UNKNOWN IMPORTED )
    set_target_properties( EOSPAC::eospac PROPERTIES
      IMPORTED_LOCATION                 "${EOSPAC_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES     "${EOSPAC_INCLUDE_DIRS}"
      IMPORTED_LINK_INTERFACE_LANGUAGES "C" )
  endif()
endif()

#------------------------------------------------------------------------------#
# End FindEOSPAC.cmake
#------------------------------------------------------------------------------#
