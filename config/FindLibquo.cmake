#.rst:
# FindLibquo
# ---------
#
# Find Libquo includes and libraries.
#
# QUO (as in "status quo") is a runtime library that aids in accommodating
# thread-level heterogeneity in dynamic, phased MPI+X applications comprising
# single- and multi-threaded libraries.
#
# https://github.com/lanl/libquo
#
# Imported Targets
# ^^^^^^^^^^^^^^^^
#
# If Libquo is found, this module defines the following :prop_tgt:`IMPORTED`
# targets::
#
#  LIBQUO::libquo         - The Libquo library.
#
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# This module will set the following variables in your project::
#
#  LIBQUO_FOUND          - True if Libquo found on the local system
#  LIBQUO_INCLUDE_DIRS   - Location of Libquo header files.
#  LIBQUO_LIBRARIES      - The Libquo libraries.
#  LIBQUO_VERSION        - The version of the discovered Libquo install.
#
# Hints
# ^^^^^
#
# Set ``LIBQUO_ROOT_DIR`` or ``QUO_HOME`` to a directory that contains a Libquo
# installation.
#
# This script expects to find libraries at ``$LIBQUO_ROOT_DIR/lib`` and the
# Libquo headers at ``$LIBQUO_ROOT_DIR/include``.  The library directory may
# optionally provide Release and Debug folders.
#
# Cache Variables
# ^^^^^^^^^^^^^^^
#
# This module may set the following variables depending on platform and type of
# Libquo installation discovered.  These variables may optionally be set to help
# this module find the correct files::
#
#  LIBQUO_LIBRARY        - Location of the Libquo library.
#  LIBQUO_LIBRARY_DEBUG  - Location of the debug Libquo library (if any).
#

#=============================================================================
# Copyright (C) 2019 Triad National Security, LLC., All rights reserved.
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

# Include these modules to handle the QUIETLY and REQUIRED arguments.
include(FindPackageHandleStandardArgs)

#=============================================================================
# If the user has provided ``LIBQUO_ROOT_DIR``, use it!  Choose items found
# at this location over system locations.
if( EXISTS "$ENV{LIBQUO_ROOT_DIR}" )
  file( TO_CMAKE_PATH "$ENV{LIBQUO_ROOT_DIR}" LIBQUO_ROOT_DIR )
elseif( EXISTS "$ENV{QUO_HOME}" )
  file( TO_CMAKE_PATH "$ENV{QUO_HOME}" LIBQUO_ROOT_DIR )
endif()
if( LIBQUO_ROOT_DIR )
  set( LIBQUO_ROOT_DIR "${LIBQUO_ROOT_DIR}" CACHE PATH
    "Prefix for Libquo installation." )
endif()

#=============================================================================
# Set LIBQUO_INCLUDE_DIRS and LIBQUO_LIBRARIES. Try to find the libraries at
# $LIBQUO_ROOT_DIR (if provided) or in standard system locations.  These
# find_library and find_path calls will prefer custom locations over standard
# locations (HINTS).  If the requested file is not found at the HINTS location,
# standard system locations will be still be searched (/usr/lib64 (Redhat),
# lib/i386-linux-gnu (Debian)).

find_path( LIBQUO_INCLUDE_DIR
  NAMES quo.h
  HINTS ${LIBQUO_ROOT_DIR}/include ${LIBQUO_INCLUDEDIR}
  PATH_SUFFIXES Release Debug
)
find_library( LIBQUO_LIBRARY
  NAMES quo
  HINTS ${LIBQUO_ROOT_DIR}/lib ${LIBQUO_LIBDIR}
  PATH_SUFFIXES Release Debug
)
# Do we also have debug versions?
find_library( LIBQUO_LIBRARY_DEBUG
  NAMES quo
  HINTS ${LIBQUO_ROOT_DIR}/lib ${LIBQUO_LIBDIR}
  PATH_SUFFIXES Debug
)
set( LIBQUO_INCLUDE_DIRS ${LIBQUO_INCLUDE_DIR} )
set( LIBQUO_LIBRARIES ${LIBQUO_LIBRARY} )

# Try to find the version.
if( NOT LIBQUO_VERSION )
  if( EXISTS "${LIBQUO_LIBDIR}/pkgconfig/libquo.pc" )
    file( STRINGS "${LIBQUO_LIBDIR}/pkgconfig/libquo.pc" libquo_version
        REGEX "Version:" )
    string( REGEX REPLACE ".*([0-9]+).([0-9]+)" "\\1" LIBQUO_MAJOR ${libquo_version} )
    string( REGEX REPLACE ".*([0-9]+).([0-9]+)" "\\2" LIBQUO_MINOR ${libquo_version} )
#    string( REGEX REPLACE ".*([0-9]+)" "\\1" LIBQUO_SUBMINOR ${libquo_h_subminor} )
  endif()
  # We might also try scraping the directory name for a regex match
  # "libquo-X.X.X"
endif()

#=============================================================================
# handle the QUIETLY and REQUIRED arguments and set LIBQUO_FOUND to TRUE if
# all listed variables are TRUE.
find_package_handle_standard_args( Libquo
  FOUND_VAR
    LIBQUO_FOUND
  REQUIRED_VARS
    LIBQUO_INCLUDE_DIR
    LIBQUO_LIBRARY
#  VERSION_VAR
#    LIBQUO_VERSION
    )

mark_as_advanced( LIBQUO_ROOT_DIR LIBQUO_VERSION LIBQUO_LIBRARY LIBQUO_INCLUDE_DIR
  LIBQUO_LIBRARY_DEBUG LIBQUO_USE_PKGCONFIG LIBQUO_CONFIG )

#=============================================================================
# Register imported libraries:
# 1. If we can find a Windows .dll file (or if we can find both Debug and
#    Release libraries), we will set appropriate target properties for these.
# 2. However, for most systems, we will only register the import location and
#    include directory.

# Look for dlls, or Release and Debug libraries.
if(WIN32)
  string( REPLACE ".lib" ".dll" LIBQUO_LIBRARY_DLL
    "${LIBQUO_LIBRARY}" )
  string( REPLACE ".lib" ".dll" LIBQUO_LIBRARY_DEBUG_DLL
    "${LIBQUO_LIBRARY_DEBUG}" )
endif()

if( LIBQUO_FOUND AND NOT TARGET LIBQUO::libquo )
  if( EXISTS "${LIBQUO_LIBRARY_DLL}" )

    # Windows systems with dll libraries.
    add_library( LIBQUO::libquo SHARED IMPORTED )

    # Windows with dlls, but only Release libraries.
    set_target_properties( LIBQUO::libquo PROPERTIES
      IMPORTED_LOCATION_RELEASE         "${LIBQUO_LIBRARY_DLL}"
      IMPORTED_IMPLIB                   "${LIBQUO_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES     "${LIBQUO_INCLUDE_DIRS}"
      IMPORTED_CONFIGURATIONS           Release
      IMPORTED_LINK_INTERFACE_LANGUAGES "C" )

    # If we have both Debug and Release libraries
    if( EXISTS "${LIBQUO_LIBRARY_DEBUG_DLL}" )
      set_property( TARGET LIBQUO::libquo APPEND PROPERTY
        IMPORTED_CONFIGURATIONS Debug )
      set_target_properties( LIBQUO::libquo PROPERTIES
        IMPORTED_LOCATION_DEBUG           "${LIBQUO_LIBRARY_DEBUG_DLL}"
        IMPORTED_IMPLIB_DEBUG             "${LIBQUO_LIBRARY_DEBUG}" )
    endif()

  else()

    # For all other environments (ones without dll libraries), create the
    # imported library targets.
    add_library( LIBQUO::libquo    UNKNOWN IMPORTED )
    set_target_properties( LIBQUO::libquo PROPERTIES
      IMPORTED_LOCATION                 "${LIBQUO_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES     "${LIBQUO_INCLUDE_DIRS}"
      IMPORTED_LINK_INTERFACE_LANGUAGES "C" )
  endif()
endif()

#------------------------------------------------------------------------------#
# End pkg_config/FindLibquo.cmake
#------------------------------------------------------------------------------#
