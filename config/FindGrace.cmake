#-----------------------------*-cmake-*----------------------------------------#
# file   config/FindGrace.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   2011 September 29
# brief  Instructions for discovering the Grace vendor libraries.
# note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

# FindGrace
# ---------
#
# Find the native Grace binary, includes and library.
#
# Grace is a WYSIWYG 2D plotting tool for the X Window System and M*tif. Grace
# runs on practically any version of Unix-like OS. As well, it has been
# successfully ported to VMS, OS/2, and Win9*/NT/2000/XP (some minor
# functionality may be missing, though). Grace is a descendant of ACE/gr, also
# known as Xmgr.
# http://plasma-gate.weizmann.ac.il/Grace/
#
# Imported Targets
# ^^^^^^^^^^^^^^^^
#
# If Grace is found, this module defines the following :prop_tgt:`IMPORTED`
# targets::
#
# Grace::grace          - The grace library.
#
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# This module will set the following variables in your project::
#
#  Grace_FOUND          - True if Grace found.
#  Grace_INCLUDE_DIRS   - where to find Grace.h, etc.
#  Grace_LIBRARIES      - List of libraries when using Grace.
#  Grace_EXECUTABLE     - main program
#
# Hints
# ^^^^^
#
# Set ``GRACE_HOME`` to a directory that contains a Grace installation.
#
# This script expects to find libraries at ``$GRACE_HOME/lib`` and the Grace
# headers at ``$GRACE_HOME/include``. The library directory may optionally
# provide Release and Debug folders.
#
# Cache Variables
# ^^^^^^^^^^^^^^^
#
# This module may set the following variables depending on platform and type of
# Grace installation discovered.  These variables may optionally be set to help
# this module find the correct files::
#
#  Grace_EXECUTABLE     - Full path to the Grace program.
#  Grace_LIBRARY        - Location of the Grace library.
#  Grace_LIBRARY_DEBUG  - Location of the debug Grace library (if any).
#  Grace_INCLUDE_DIRS   - Location of the Grace headers.
#

#=============================================================================
# Include these modules to handle the QUIETLY and REQUIRED arguments.
include(FindPackageHandleStandardArgs)

#=============================================================================
# If the user has provided ``GRACE_HOME``, use it!  Choose items found at
# this location over system locations.
if( EXISTS "$ENV{GRACE_HOME}" )
  file( TO_CMAKE_PATH "$ENV{GRACE_HOME}" GRACE_HOME )
  set( GRACE_HOME "${GRACE_HOME}" CACHE PATH
    "Prefix for Grace installation." )
endif()

# If environment variables are not set (GRACE_INC_DIR), look for the binary and
# try to guess appropriate location of library and headers.
find_program( Grace_EXECUTABLE NAMES xmgrace ggrace)

if( EXISTS ${Grace_EXECUTABLE} AND NOT IS_DIRECTORY "${GRACE_HOME}" )
  get_filename_component( GRACE_HOME ${Grace_EXECUTABLE} DIRECTORY )
  get_filename_component( GRACE_HOME "${GRACE_HOME}/.." REALPATH )
endif()

#=============================================================================
# Set Grace_INCLUDE_DIRS and Grace_LIBRARIES. Try to find the libraries at
# $GRACE_HOME (if provided) or in standard system locations.  These find_library
# and find_path calls will prefer custom locations over standard locations
# (HINTS).  If the requested file is not found at the HINTS location, standard
# system locations will be still be searched (/usr/lib64 (Redhat),
# lib/i386-linux-gnu (Debian)).

find_path( Grace_INCLUDE_DIR
  NAMES grace_np.h
  HINTS ${GRACE_HOME}/include
  PATH_SUFFIXES Release Debug
  )

if( WIN32 )
  set( Grace_LIBRARY_NAME libgrace_dll.lib)
elseif( APPLE )
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".dylib")
else()
  set( Grace_LIBRARY_NAME grace_np)
endif()

find_library(Grace_LIBRARY
  NAMES ${Grace_LIBRARY_NAME}
  HINTS ${GRACE_HOME}/lib
  PATH_SUFFIXES Release Debug
  )
find_library(Grace_LIBRARY_DEBUG
  NAMES ${Grace_LIBRARY_NAME}
  HINTS ${GRACE_HOME}/lib
  PATH_SUFFIXES Debug
  )

set( Grace_INCLUDE_DIRS ${Grace_INCLUDE_DIR} )
set( Grace_LIBRARIES ${Grace_LIBRARY} )

# Try to find the version.
if( NOT Grace_VERSION )
  if( EXISTS "${Grace_EXECUTABLE}" )
    execute_process(
      COMMAND ${Grace_EXECUTABLE} -version
      OUTPUT_VARIABLE grace_ver_out
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_STRIP_TRAILING_WHITESPACE
      )
    string( REPLACE "\n" ";"  grace_ver_out "${grace_ver_out}" )
    foreach( line ${grace_ver_out} )
      if( "${line}" MATCHES "Grace-" )
        set( Grace_VERSION ${line} )
      endif()
    endforeach()
    string( REGEX REPLACE ".*([0-9]+)[.]([0-9]+)[.]([0-9]+).*" "\\1"
      Grace_VER_MAJOR "${Grace_VERSION}" )
    string( REGEX REPLACE ".*([0-9]+)[.]([0-9]+)[.]([0-9]+).*" "\\2"
      Grace_VER_MINOR "${Grace_VERSION}" )
  endif()
endif()

#=============================================================================
# handle the QUIETLY and REQUIRED arg
find_package_handle_standard_args(Grace
  FOUND_VAR     Grace_FOUND
  REQUIRED_VARS Grace_INCLUDE_DIR Grace_LIBRARY
  VERSION_VAR   Grace_VERSION
  )

mark_as_advanced( GRACE_HOME Grace_VERSION Grace_LIBRARY Grace_LIBRARY_DEBUG
  Grace_INCLUDE_DIR Grace_EXECUTABLE)

#=============================================================================
# Register imported libraries:
# 1. If we can find a Windows .dll file (or if we can find both Debug and
#    Release libraries), we will set appropriate target properties for these.
# 2. However, for most systems, we will only register the import location and
#    include directory.

# Look for dlls, or Release and Debug libraries.
if(WIN32)
  string( REPLACE ".lib" ".dll" Grace_LIBRARY_DLL
    "${Grace_LIBRARY}" )
  string( REPLACE ".lib" ".dll" Grace_LIBRARY_DEBUG_DLL
    "${Grace_LIBRARY_DEBUG}" )
endif()

if( Grace_FOUND AND NOT TARGET Grace::libgrace )
  if( EXISTS "${Grace_LIBRARY_DLL}" )

    # Windows systems with dll libraries.
    add_library( Grace::libgrace SHARED IMPORTED )

    # Windows with dlls, but only Release libraries.
    set_target_properties( Grace::libgrace PROPERTIES
      IMPORTED_LOCATION_RELEASE         "${Grace_LIBRARY_DLL}"
      IMPORTED_IMPLIB                   "${Grace_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES     "${Grace_INCLUDE_DIRS}"
      IMPORTED_CONFIGURATIONS           Release
      IMPORTED_LINK_INTERFACE_LANGUAGES "C" )

    # If we have both Debug and Release libraries
    if( EXISTS "${Grace_LIBRARY_DEBUG_DLL}" )
      set_property( TARGET Grace::libgrace APPEND PROPERTY
        IMPORTED_CONFIGURATIONS Debug )
      set_target_properties( Grace::libgrace PROPERTIES
        IMPORTED_LOCATION_DEBUG           "${Grace_LIBRARY_DEBUG_DLL}"
        IMPORTED_IMPLIB_DEBUG             "${Grace_LIBRARY_DEBUG}" )
    endif()

  else()

    # For all other environments (ones without dll libraries), create the
    # imported library targets.
    add_library( Grace::libgrace    UNKNOWN IMPORTED )
    set_target_properties( Grace::libgrace PROPERTIES
      IMPORTED_LOCATION                 "${Grace_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES     "${Grace_INCLUDE_DIRS}"
      IMPORTED_LINK_INTERFACE_LANGUAGES "C" )
  endif()
endif()

#------------------------------------------------------------------------------#
# End of FindGrace.cmake
#------------------------------------------------------------------------------#
