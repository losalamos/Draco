#.rst:
# FindSuperLU_DIST
# --------
#
# Find the SuperLU_DIST includes and libraries.
#
# SuperLU is a general purpose library for the direct solution of large, sparse,
# nonsymmetric systems of linear equations on high performance machines. The
# library is written in C and is callable from either C or Fortran. The library
# routines will perform an LU decomposition with partial pivoting and triangular
# system solves through forward and back substitution. The LU factorization
# routines can handle non-square matrices but the triangular solves are
# performed only for square matrices. The matrix columns may be preordered
# (before factorization) either through library or user supplied routines. This
# preordering for sparsity is completely separate from the
# factorization. Working precision iterative refinement subroutines are provided
# for improved backward stability. Routines are also provided to equilibrate the
# system, estimate the condition number, calculate the relative backward error,
# and estimate error bounds for the refined solutions.
# http://crd-legacy.lbl.gov/~xiaoye/SuperLU/
#
# Imported Targets
# ^^^^^^^^^^^^^^^^
#
# If SuperLU_DIST is found, this module defines the following
# :prop_tgt:`IMPORTED` targets::
#
# SuperLU_DIST::superludist      - The main SuperLU_DIST library.
#
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# This module will set the following variables in your project::
#
#  SuperLU_DIST_FOUND          - True if SuperLU_DIST found on the local system
#  SuperLU_DIST_INCLUDE_DIRS   - Location of SuperLU_DIST header files.
#  SuperLU_DIST_LIBRARIES      - The SuperLU_DIST libraries.
#  SuperLU_DIST_VERSION        - The version of the discovered SuperLU_DIST
#                                install.
#
# Hints
# ^^^^^
#
# Set ``SuperLU_DIST_ROOT_DIR`` to a directory that contains a SuperLU_DIST
# installation.
#
# This script expects to find libraries at ``$SuperLU_DIST_ROOT_DIR/lib`` and
# the SuperLU_DIST headers at ``$SuperLU_DIST_ROOT_DIR/include``.  The library
# directory may optionally provide Release and Debug folders.
#
# Cache Variables
# ^^^^^^^^^^^^^^^
#
# This module may set the following variables depending on platform and type of
# SuperLU_DIST installation discovered.  These variables may optionally be set
# to help this module find the correct files::
#
#  SuperLU_DIST_LIBRARY       - Location of the SuperLU_DIST library.
#  SuperLU_DIST_LIBRARY_DEBUG - Location of the debug SuperLU_DIST library
#                               (if any).
#

#=============================================================================
# Copyright 2016 Kelly Thompson <kgt@lanl.gov>
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
# If the user has provided ``SuperLU_DIST_ROOT_DIR``, use it!  Choose items
# found at this location over system locations.
if( EXISTS "$ENV{SuperLU_DIST_ROOT_DIR}" )
  file( TO_CMAKE_PATH "$ENV{SuperLU_DIST_ROOT_DIR}" SuperLU_DIST_ROOT_DIR )
  set( SuperLU_DIST_ROOT_DIR "${SuperLU_DIST_ROOT_DIR}" CACHE PATH
    "Prefix for SuperLU_DIST installation." )
endif()

#=============================================================================
# Set SuperLU_DIST_INCLUDE_DIRS and SuperLU_DIST_LIBRARIES. Try to find the
# libraries at $SuperLU_DIST_ROOT_DIR (if provided) or in standard system
# locations.  These find_library and find_path calls will prefer custom
# locations over standard locations (HINTS).  If the requested file is not found
# at the HINTS location, standard system locations will be still be searched
# (/usr/lib64 (Redhat), lib/i386-linux-gnu (Debian)).

find_path( SuperLU_DIST_INCLUDE_DIR
  NAMES superlu_defs.h
  HINTS ${SuperLU_DIST_ROOT_DIR}/include ${SuperLU_DIST_INCLUDEDIR}
  PATH_SUFFIXES Release Debug
  )
find_library( SuperLU_DIST_LIBRARY
  NAMES superlu_dist superludist
  HINTS ${SuperLU_DIST_ROOT_DIR}/lib ${SuperLU_DIST_LIBDIR}
  PATH_SUFFIXES Release Debug
  )
# Do we also have debug versions?
find_library( SuperLU_DIST_LIBRARY_DEBUG
  NAMES superlu_dist superludist
  HINTS ${SuperLU_DIST_ROOT_DIR}/lib ${SuperLU_DIST_LIBDIR}
  PATH_SUFFIXES Debug
  )
set( SuperLU_DIST_INCLUDE_DIRS ${SuperLU_DIST_INCLUDE_DIR} )
set( SuperLU_DIST_LIBRARIES ${SuperLU_DIST_LIBRARY} )

# Try to find the version.
# Not sure if there is an easy way to query headers or libs for the version tag.
if( NOT SuperLU_DIST_VERSION )
  # Try scraping the directory name for a regex match "SuperLU_DIST-X.X.X"
  string( REGEX REPLACE ".*SuperLU_DIST-([0-9][.][0-9][.][0-9])" "\\1"
    SuperLU_DIST_MAJOR ${SuperLU_DIST_INCLUDE_DIR} )
  string( REGEX REPLACE ".*SuperLU_DIST-([0-9][.][0-9][.][0-9])" "\\2"
    SuperLU_DIST_MAJOR ${SuperLU_DIST_INCLUDE_DIR} )
  string( REGEX REPLACE ".*SuperLU_DIST-([0-9][.][0-9][.][0-9])" "\\3"
    SuperLU_DIST_MAJOR ${SuperLU_DIST_INCLUDE_DIR} )
  set( SuperLU_DIST_VERSION "${SuperLU_DIST_MAJOR}.${SuperLU_DIST_MINOR}.${SuperLU_DIST_SUBMINOR}" )
  if( ${SuperLU_DIST_VERSION} STREQUAL ".." )
    set( SuperLU_DIST_VERSION "0.0.0")
  endif()
endif()

#=============================================================================
# handle the QUIETLY and REQUIRED arguments and set SuperLU_DIST_FOUND to TRUE if
# all listed variables are TRUE.
find_package_handle_standard_args( SuperLU_DIST
  FOUND_VAR
    SuperLU_DIST_FOUND
  REQUIRED_VARS
    SuperLU_DIST_INCLUDE_DIR
    SuperLU_DIST_LIBRARY
  VERSION_VAR
    SuperLU_DIST_VERSION
  )

mark_as_advanced( SuperLU_DIST_ROOT_DIR SuperLU_DIST_VERSION SuperLU_DIST_LIBRARY
  SuperLU_DIST_INCLUDE_DIR SuperLU_DIST_LIBRARY_DEBUG )

#=============================================================================
# Register imported libraries:
# 1. If we can find a Windows .dll file (or if we can find both Debug and
#    Release libraries), we will set appropriate target properties for these.
# 2. However, for most systems, we will only register the import location and
#    include directory.

# Look for dlls, or Release and Debug libraries.
if(WIN32)
  string( REPLACE ".lib" ".dll" SuperLU_DIST_LIBRARY_DLL
    "${SuperLU_DIST_LIBRARY}" )
  string( REPLACE ".lib" ".dll" SuperLU_DIST_LIBRARY_DEBUG_DLL
    "${SuperLU_DIST_LIBRARY_DEBUG}" )
endif()

if( SuperLU_DIST_FOUND AND NOT TARGET SuperLU_DIST::superludist )
  if( EXISTS "${SuperLU_DIST_LIBRARY_DLL}" )

    # Windows systems with dll libraries.
    add_library( SuperLU_DIST::superludist SHARED IMPORTED )

    # Windows with dlls, but only Release libraries.
    set_target_properties( SuperLU_DIST::superludist PROPERTIES
      IMPORTED_LOCATION_RELEASE         "${SuperLU_DIST_LIBRARY_DLL}"
      IMPORTED_IMPLIB                   "${SuperLU_DIST_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES     "${SuperLU_DIST_INCLUDE_DIRS}"
      IMPORTED_CONFIGURATIONS           Release
      IMPORTED_LINK_INTERFACE_LANGUAGES "C" )

    # If we have both Debug and Release libraries
    if( EXISTS "${SuperLU_DIST_LIBRARY_DEBUG_DLL}" )
      set_property( TARGET SuperLU_DIST::superludist APPEND PROPERTY
        IMPORTED_CONFIGURATIONS Debug )
      set_target_properties( SuperLU_DIST::superludist PROPERTIES
        IMPORTED_LOCATION_DEBUG "${SuperLU_DIST_LIBRARY_DEBUG_DLL}"
        IMPORTED_IMPLIB_DEBUG   "${SuperLU_DIST_LIBRARY_DEBUG}" )
    endif()

  else()

    # For all other environments (ones without dll libraries), create the
    # imported library targets.
    add_library( SuperLU_DIST::superludist UNKNOWN IMPORTED )
    set_target_properties( SuperLU_DIST::superludist PROPERTIES
      IMPORTED_LOCATION                 "${SuperLU_DIST_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES     "${SuperLU_DIST_INCLUDE_DIRS}"
      IMPORTED_LINK_INTERFACE_LANGUAGES "C" )
  endif()
endif()
