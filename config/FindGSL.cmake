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
#  GSL_INCLUDE_DIR    - Locaiton of GSL header files.
#  GSL_LIBRARIES      - List of GSL libraries.
#  GSL_VERSION        - 
#  GSL_ROOT_DIR       - 
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
# ::
#
#    1. Set GSL_ROOT_DIR to a directory that contains bin/gsl-config.
#    2. Set GSL_INC_DIR and GSL_LIB_DIR to the locations of the
#       installed GSL include files and libraries are located.
#
# When configuration is successful, GSL_INCLUDE_DIR will be set to
# the location of the GSL header files and GSL_LIBRARIES will be set
# to a list of fully qualified path locations of the GSL libraries. 

#=============================================================================
# Copyright (C) 2010-2014 Los Alamos National Security, LLC.
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================

# include this to handle the QUIETLY and REQUIRED arguments
include(FindPackageHandleStandardArgs)
include(GetPrerequisites)

# In the future, we might be able to modify this module to use:
# find_package(PkgConfig)
# pkg_check_modules( PC_GSL QUIET GSL )

#=============================================================================

# If the user has provided $GSL_ROOT_DIR, use it!  Choose items found
# at this location over system locations.
if( EXISTS "$ENV{GSL_ROOT_DIR}" )
  set( GSL_ROOT_DIR "$ENV{GSL_ROOT_DIR}" CACHE PATH "Prefix for GSL installation." )
endif()

# Find and use 'gsl-config' to determine version, library locations, etc.
find_program( GSL_CONFIG gsl-config
  HINTS ${GSL_ROOT_DIR}/bin
  )

if( EXISTS "${GSL_CONFIG}" )

  # If we have gsl-config, use it to locate the include directory, libraries 
  # and version number.  Some installations of GSL (including Win32) will not 
  # have this script so we must rely on the libraries being in system locations 
  # or the user providing GSL_BOOST_ROOT.
  
  # set the version:
  exec_program( "${GSL_CONFIG}"
    ARGS --version
    OUTPUT_VARIABLE GSL_VERSION )

  # set CXXFLAGS to be fed into CXX_FLAGS by the user:
  exec_program( "${GSL_CONFIG}"
    ARGS --cflags
    OUTPUT_VARIABLE GSL_CFLAGS )
  
  # set INCLUDE_DIR to prefix+include
  exec_program( "${GSL_CONFIG}"
    ARGS --prefix
    OUTPUT_VARIABLE GSL_PREFIX)
  set( GSL_INCLUDE_DIR ${GSL_PREFIX}/include CACHE STRING INTERNAL )
  # Reset GSL_ROOT_DIR based on what we found.
  set( GSL_ROOT_DIR ${GSL_PREFIX} CACHE STRING "Location of GSL.")
  
  # set link libraries and link flags
  exec_program( "${GSL_CONFIG}"
    ARGS --libs
    OUTPUT_VARIABLE GSL_LIBRARIES )

  # extract link dirs (if any)
  # use regular expression to match wildcard equivalent "-L*<endchar>" with 
  # <endchar> is a space or a semicolon
  string( REGEX MATCHALL "[-][L]([^ ;])+" GSL_LINK_DIRECTORIES_WITH_PREFIX 
    "${GSL_LIBRARIES}" )
  if( GSL_LINK_DIRECTORIES_WITH_PREFIX )
    string( REGEX REPLACE "[-][L]" "" GSL_LINK_DIRECTORIES 
      ${GSL_LINK_DIRECTORIES_WITH_PREFIX} )
  endif()

endif()

#=============================================================================

if( NOT EXISTS "${GSL_INCLUDE_DIR}" )
  find_path( GSL_INCLUDE_DIR 
    NAMES gsl/gsl_sf.h
    HINTS ${GSL_ROOT_DIR}/include
    )
endif()

find_library( GSL_LIBRARY
  NAMES gsl 
  HINTS 
  ${GSL_ROOT_DIR}/lib
  ${GSL_LINK_DIRECTORIES}
  PATH_SUFFIXES Release Debug
  )

# On Linux systems, this library is typically named libgslcblas.a
# On Windows systems, this library is typically named cblas.dll
find_library( GSL_CBLAS_LIBRARY
  NAMES gslcblas cblas
  HINTS
  ${GSL_ROOT_DIR}/lib
  ${GSL_LINK_DIRECTORIES}
  PATH_SUFFIXES Release Debug
  )

# For windows platforms, look for the debug libraries also
if(WIN32)
  find_library( gsl_library_debug
    names gsl 
    hints ${gsl_root_dir}/lib
    path_suffixes debug
    )
  find_library( gsl_cblas_library_debug
    names gslcblas cblas
    hints ${gsl_root_dir}/lib
    path_suffixes debug
    )
endif()

# handle the QUIETLY and REQUIRED arguments and set GSL_FOUND to TRUE if 
# all listed variables are TRUE
find_package_handle_standard_args( GSL
  FOUND_VAR GSL_FOUND
  REQUIRED_VARS 
  GSL_INCLUDE_DIR 
  GSL_LIBRARY
  GSL_CBLAS_LIBRARY 
  GSL_ROOT_DIR
  VERSION_VAR GSL_VERSION
  )
mark_as_advanced( GSL_FOUND GSL_LIBRARY GSL_CBLAS_LIBRARY
  GSL_INCLUDE_DIR GSL_ROOT_DIR )

# Register imported libraries
if( GSL_FOUND )

  # If this platform doesn't support shared libraries (e.g. cross
  # compiling), assume static. This suppresses cmake (3.0.0) warnings
  # of the form: 
  #     "ADD_LIBRARY called with MODULE option but the target platform
  #     does not support dynamic linking.  Building a STATIC library
  #     instead.  This may lead to problems."
  if( TARGET_SUPPORTS_SHARED_LIBS )
    add_library( gsl::gsl      UNKNOWN IMPORTED )
    add_library( gsl::gslcblas UNKNOWN IMPORTED )
  else()
    add_library( gsl::gsl      STATIC IMPORTED )
    add_library( gsl::gslcblas STATIC IMPORTED )
  endif()
  if( EXISTS "${GSL_LIBRARY_DEBUG}" )
    set_target_properties( gsl::gslcblas PROPERTIES
      IMPORTED_LOCATION_RELEASE "${GSL_CBLAS_LIBRARY}"
      IMPORTED_LOCATION_DEBUG "${GSL_CBLAS_LIBRARY_DEBUG}"
      INTERFACE_INCLUDE_DIRECTORIES "${GSL_INCLUDE_DIR}"
      IMPORTED_LINK_INTERFACE_LANGUAGES "C" )
    set_target_properties( gsl::gsl PROPERTIES
      IMPORTED_LOCATION_RELEASE "${GSL_LIBRARY}"
      IMPORTED_LOCATION_DEBUG "${GSL_LIBRARY_DEBUG}"
      INTERFACE_INCLUDE_DIRECTORIES "${GSL_INCLUDE_DIR}" 
      IMPORTED_LINK_INTERFACE_LANGUAGES "C"
      IMPORTED_LINK_INTERFACE_LIBRARIES gsl::gslcblas )        
    # cleanup
    unset( GSL_LIBRARY_DEBUG CACHE)
    unset( GSL_CBLAS_LIBRARY_DEBUG CACHE)
  else()
    set_target_properties( gsl::gslcblas PROPERTIES
      IMPORTED_LOCATION "${GSL_CBLAS_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${GSL_INCLUDE_DIR}"
      IMPORTED_LINK_INTERFACE_LANGUAGES "C" )
    set_target_properties( gsl::gsl PROPERTIES
      IMPORTED_LOCATION "${GSL_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${GSL_INCLUDE_DIR}" 
      IMPORTED_LINK_INTERFACE_LANGUAGES "C"
      IMPORTED_LINK_INTERFACE_LIBRARIES gsl::gslcblas )
  endif()

endif()

# Include some information that can be printed by the build system.
include( FeatureSummary )
set_package_properties( GSL PROPERTIES
  DESCRIPTION "Gnu Scientific Library"
  URL "www.gnu.org/software/gsl"
  PURPOSE "The GNU Scientific Library (GSL) is a numerical library for C and C++ programmers."
  )

# cleanup
unset( GSL_CONFIG CACHE )

# See also
# http://sourceforge.net/p/openmodeller/svn/HEAD/tree/trunk/openmodeller/cmake/FindGSL.cmake
# http://www.cmake.org/cmake/help/git-master/manual/cmake-developer.7.html#modules

