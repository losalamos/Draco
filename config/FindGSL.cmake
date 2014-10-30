#.rst
# Find GSL
# --------
#
# Find the native GSL includes and library
#
# The GNU Scientific Library (GSL) is a numerical library for C and
# C++ programmers. It is free software under the GNU General Public
# License. 
#
# === Variables ===
#
# This module will set the following variables per language in your
# project:
#
# ::
#
#  GSL_FOUND          - True if GSL found on the local system
#  GSL_INCLUDE_DIR    - Location of GSL header files.
#  GSL_LIBRARY        - The GSL library.
#  GSL_CBLAS_LIBRARY  - The GSL CBLAS library.
#  GSL_VERSION        - The version of the discovered GSL install.
#  GSL_ROOT_DIR       - The top level directory of the discovered GSL 
#                       install (useful if GSL is not in a standard location)
#
# It will also provide the cmake library target names gsl::gsl and
# gsl::gslcblas. 
#
# === Usage ===
#
# To use this module, simply call FindGSL from a CMakeLists.txt file, or
# run find_package(GSL), then run CMake. If you are happy with the
# auto- detected configuration for your language, then you're done.  If
# not, you have one option:
#
#    - Set GSL_ROOT_DIR to a directory that contains a GSL install. This 
#      script expects to find libraries at $GSL_ROOT_DIR/lib and the 
#      GSL headers at $GSL_ROOT_DIR/include/gsl.  The library directory
#      may optionally provide Release and Debug folders.  For Unix-like 
#      systems, this script will use $GSL_ROOT_DIR/bin/gsl-config (if found)
#      to aid in the discovery GSL.
#
# When configuration is successful, GSL_INCLUDE_DIR will be set to
# the location of the GSL header files and GSL_LIBRARY and GSL_CBLAS_LIBRARY 
# will be set to fully qualified path locations of the GSL libraries. 

#=============================================================================
# Copyright (C) 2010-2014 Los Alamos National Security, LLC.
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================

# Warn about using with older versions of cmake...
# - This module uses PkgConfig (since 2.8), import libraries (since
#   2.6), but developed and only tested with cmake 3.0+
if(CMAKE_MINIMUM_REQUIRED_VERSION VERSION_LESS 3.0.0)
  message(AUTHOR_WARNING
    "Your project should require at least CMake 3.0.0 to use FindFoo.cmake")
endif()

# Include these modules to handle the QUIETLY and REQUIRED arguments. See
# - http://www.cmake.org/cmake/help/v3.0/module/FindPackageHandleStandardArgs.html
# - http://www.cmake.org/cmake/help/v3.0/module/GetPrerequisites.html
include(FindPackageHandleStandardArgs)
include(GetPrerequisites)

#=============================================================================
# If the user has provided $GSL_ROOT_DIR, use it!  Choose items found
# at this location over system locations.
if( EXISTS "$ENV{GSL_ROOT_DIR}" )
  set( GSL_ROOT_DIR "$ENV{GSL_ROOT_DIR}" CACHE PATH "Prefix for GSL installation." )
endif()
if( NOT EXISTS "${GSL_ROOT_DIR}" )
  set( GSL_USE_PKGCONFIG ON )
endif()

#=============================================================================
# As a first try, use the PkgConfig module.  This will work on many
# *NIX systems.  See http://www.cmake.org/cmake/help/v3.0/module/FindPkgConfig.html
if( GSL_USE_PKGCONFIG )
  find_package(PkgConfig)
  pkg_check_modules( GSL QUIET gsl )

  if( EXISTS "${GSL_INCLUDEDIR}" )
    get_filename_component( GSL_ROOT_DIR "${GSL_INCLUDEDIR}" DIRECTORY CACHE)
  endif() 
endif()

#=============================================================================
# Set GSL_INCLUDE_DIR, GSL_LIBRARY and GSL_CBLAS_LIBRARY.
# If we skipped the PkgConfig step, try to find the libraries at
# $GSL_ROOT_DIR (if provided) or in standard system locations.
  
find_path( GSL_INCLUDE_DIR 
  NAMES gsl/gsl_sf.h
  HINTS ${GSL_ROOT_DIR}/include ${GSL_INCLUDEDIR}
)
find_library( GSL_LIBRARY
  NAMES gsl 
  HINTS ${GSL_ROOT_DIR}/lib ${GSL_LIBDIR}
  PATH_SUFFIXES Release Debug
)
find_library( GSL_CBLAS_LIBRARY
  NAMES gslcblas cblas
  HINTS ${GSL_ROOT_DIR}/lib ${GSL_LIBDIR}
  PATH_SUFFIXES Release Debug
)
# Do we also have debug versions?
find_library( GSL_LIBRARY_DEBUG
  NAMES gsl 
  HINTS ${GSL_ROOT_DIR}/lib ${GSL_LIBDIR}
  PATH_SUFFIXES Debug
)
find_library( GSL_CBLAS_LIBRARY_DEBUG
  NAMES gslcblas cblas
  HINTS ${GSL_ROOT_DIR}/lib ${GSL_LIBDIR}
  PATH_SUFFIXES Debug
)

# If we didn't use PkgConfig, try to find the version via gsl-config
# or by reading gsl_version.h.
if( NOT GSL_VERSION )
  # 1. If gsl-config exists, query for the version.
  find_program( GSL_CONFIG
    NAMES gsl-config
    HINTS "${GSL_ROOT_DIR}/bin"
    )
  if( EXISTS "${GSL_CONFIG}" )
    exec_program( "${GSL_CONFIG}"
      ARGS --version
      OUTPUT_VARIABLE GSL_VERSION )
  endif()
  
  # 2. If gsl-config is not available, try looking in gsl/gsl_version.h
  if( NOT GSL_VERSION AND EXISTS "${GSL_INCLUDE_DIR}/gsl/gsl_version.h" )
    file( STRINGS "${GSL_INCLUDE_DIR}/gsl/gsl_version.h" gsl_version_h_contents )
    foreach( entry ${gsl_version_h_contents} )
      if( ${entry} MATCHES "define GSL_VERSION")
         string( REGEX REPLACE ".*([0-9].[0-9][0-9]).*" "\\1" GSL_VERSION ${entry} )
      endif()
    endforeach()    
  endif()
  
  # might also try scraping the directory name for a regex match "gsl-X.X"
endif()

#=============================================================================
# handle the QUIETLY and REQUIRED arguments and set GSL_FOUND to TRUE if 
# all listed variables are TRUE
find_package_handle_standard_args( GSL
  FOUND_VAR
    GSL_FOUND
  REQUIRED_VARS 
    GSL_INCLUDE_DIR 
    GSL_LIBRARY
    GSL_CBLAS_LIBRARY 
    GSL_ROOT_DIR
  VERSION_VAR
    GSL_VERSION
  )
mark_as_advanced( GSL_FOUND GSL_LIBRARY GSL_CBLAS_LIBRARY
  GSL_INCLUDE_DIR GSL_ROOT_DIR GSL_VERSION )

#=============================================================================
# Register imported libraries:
# 1. If we can find a Windows .dll file (or if we can find both Debug and 
#    Release libraries), we will set appropriate target properties for these.  
# 2. However, for most systems, we will only register the import location and 
#    include directory.

# Look for dlls, or Release and Debug libraries.

if(WIN32)
  string( REPLACE ".lib" ".dll" GSL_LIBRARY_DLL       "${GSL_LIBRARY}" )
  string( REPLACE ".lib" ".dll" GSL_CBLAS_LIBRARY_DLL "${GSL_CBLAS_LIBRARY}" )
  string( REPLACE ".lib" ".dll" GSL_LIBRARY_DEBUG_DLL "${GSL_LIBRARY_DEBUG}" )
  string( REPLACE ".lib" ".dll" GSL_CBLAS_LIBRARY_DEBUG_DLL "${GSL_CBLAS_LIBRARY_DEBUG}" )
endif()

if( GSL_FOUND )
  if( EXISTS "${GSL_LIBRARY_DLL}" AND EXISTS "${GSL_CBLAS_LIBRARY_DLL}")
    # Windows systems with dll libraries.
    add_library( gsl::gsl      SHARED IMPORTED )
    add_library( gsl::gslcblas SHARED IMPORTED )
    # If we have both Debug and Release libraries
    if( EXISTS "${GSL_LIBRARY_DEBUG_DLL}" AND EXISTS "${GSL_CBLAS_LIBRARY_DEBUG_DLL}")
      set_target_properties( gsl::gslcblas PROPERTIES
        IMPORTED_LOCATION                 "${GSL_CBLAS_LIBRARY_DLL}"
        IMPORTED_IMPLIB                   "${GSL_CBLAS_LIBRARY}"
        IMPORTED_LOCATION_DEBUG           "${GSL_CBLAS_LIBRARY_DEBUG_DLL}"
        IMPORTED_IMPLIB_DEBUG             "${GSL_CBLAS_LIBRARY_DEBUG}"
        INTERFACE_INCLUDE_DIRECTORIES     "${GSL_INCLUDE_DIR}"
        IMPORTED_LINK_INTERFACE_LANGUAGES "C" )
      set_target_properties( gsl::gsl PROPERTIES
        IMPORTED_LOCATION                 "${GSL_LIBRARY_DLL}"
        IMPORTED_IMPLIB                   "${GSL_LIBRARY}"
        IMPORTED_LOCATION_DEBUG           "${GSL_LIBRARY_DEBUG_DLL}"
        IMPORTED_IMPLIB_DEBUG             "${GSL_LIBRARY_DEBUG}"
        INTERFACE_INCLUDE_DIRECTORIES     "${GSL_INCLUDE_DIR}" 
        IMPORTED_LINK_INTERFACE_LANGUAGES "C"
        IMPORTED_LINK_INTERFACE_LIBRARIES gsl::gslcblas )
    else()
      # Windows with dlls, but only Release libraries.
      set_target_properties( gsl::gslcblas PROPERTIES
        IMPORTED_LOCATION                 "${GSL_CBLAS_LIBRARY_DLL}"
        IMPORTED_IMPLIB                   "${GSL_CBLAS_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES     "${GSL_INCLUDE_DIR}"
        IMPORTED_LINK_INTERFACE_LANGUAGES "C" )
      set_target_properties( gsl::gsl PROPERTIES
        IMPORTED_LOCATION                 "${GSL_LIBRARY_DLL}"
        IMPORTED_IMPLIB                   "${GSL_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES     "${GSL_INCLUDE_DIR}" 
        IMPORTED_LINK_INTERFACE_LANGUAGES "C"
        IMPORTED_LINK_INTERFACE_LIBRARIES gsl::gslcblas )    
    endif()

  else()
    # *NIX systems or simple installs.
    add_library( gsl::gsl      UNKNOWN IMPORTED )
    add_library( gsl::gslcblas UNKNOWN IMPORTED )
    set_target_properties( gsl::gslcblas PROPERTIES
      IMPORTED_LOCATION                 "${GSL_CBLAS_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES     "${GSL_INCLUDE_DIR}"
      IMPORTED_LINK_INTERFACE_LANGUAGES "C" )
    set_target_properties( gsl::gsl PROPERTIES
      IMPORTED_LOCATION                 "${GSL_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES     "${GSL_INCLUDE_DIR}" 
      IMPORTED_LINK_INTERFACE_LANGUAGES "C"
      IMPORTED_LINK_INTERFACE_LIBRARIES gsl::gslcblas )
  endif()
endif()

#=============================================================================
# Include some information that can be printed by the build system.
include( FeatureSummary )
set_package_properties( GSL PROPERTIES
  DESCRIPTION "Gnu Scientific Library"
  URL "www.gnu.org/software/gsl"
  PURPOSE "The GNU Scientific Library (GSL) is a numerical library for C and C++ programmers."
  )

#=============================================================================
# cleanup
unset( GSL_USE_PKGCONFIG )
unset( GSL_CONFIG CACHE )

#=============================================================================
# References:
# - http://www.cmake.org/cmake/help/git-master/manual/cmake-developer.7.html#modules
