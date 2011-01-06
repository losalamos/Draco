# - Find the DRACO binary suite (based on FindImageMagick.cmake)
#
# This module will search for a set of DRACO tools specified as
# components in the FIND_PACKAGE call. Typical components include, but
# are not limited to (future versions of DRACO might have additional
# components not listed here):
#
# RTT_Format_Reader
# c4
# cdi
# cdi_analytic
# cdi_eospac
# cdi_gandolf
# ds++
# fit
# fpe_trap
# lapack_wrap
# linear
# meshReaders
# mesh_element
# min
# norms
# ode
# parser
# pcgWrap
# plot2D
# quadrature
# rng
# roots
# shared_lib
# special_functions
# timestep
# traits
# units
# viz
# xm
#
# If no component is specified in the FIND_PACKAGE call, then it only
# searches for the DRACO executable directory. This code defines the
# following variables:
#
#  DRACO_FOUND                    - TRUE if all components are found.
#  DRACO_DIR                      - Full path to the DRACO install directory.
#  DRACO_LIBRARY_DIR              - Full path to libraries directory.
#  DRACO_INCLUDE_DIR              - Full path to include directory.
#  DRACO_LIBRARIES                - Full paths to all libraries.
#  DRACO_<component>_FOUND        - TRUE if <component> is found.
#  DRACO_<component>_INCLUDE_DIR  - Full path to <component> include dirs.
#  DRACO_<component>_LIBRARIES    - Full path to <component> libraries.
#
# Example Usages:
#  FIND_PACKAGE(DRACO)
#  FIND_PACKAGE(DRACO COMPONENTS ds++)
#  FIND_PACKAGE(DRACO COMPONENTS ds++ c4 REQUIRED)
#  FIND_PACKAGE(DRACO COMPONENTS ds++ QUIET)
#  FIND_PACKAGE(DRACO COMPONENTS RTT_Format_Reader meshReaders mesh_element ds++)
#
# Note that the standard FIND_PACKAGE features are supported
# (i.e., QUIET, REQUIRED, etc.).
#
# You may need to set environment or cmake variables to help the build
# system find draco
# - DRACO_DIR
# - DRACO_INCLUDE_DIR
# - DRACO_LIBRARY_DIR
#

#---------------------------------------------------------------------
# Helper functions
# e.g.:     FIND_DRACO_API(ds++ SP.hh ds++ CORE_RL_ds++_ )
#---------------------------------------------------------------------
function(FIND_DRACO_API component header)
  set(DRACO_${component}_FOUND FALSE PARENT_SCOPE)

  find_path(DRACO_${component}_INCLUDE_DIR
    NAMES Release.hh # ${header}
    PATHS
      ${DRACO_INCLUDE_DIR}
      ${DRACO_DIR}/include
      ${CMAKE_INSTALL_PREFIX}/include
#      "[HKEY_LOCAL_MACHINE\\SOFTWARE\\DRACO\\Current;BinPath]/include"
    PATH_SUFFIXES
       ${component}
#      DRACO
    DOC "Path to the DRACO <componenet> include dir."
    )
  find_library(DRACO_${component}_LIBRARY
    NAMES rtt_${component}
    PATHS
      ${DRACO_LIBRARY_DIR}
      ${DRACO_DIR}/lib
      ${CMAKE_INSTALL_PREFIX}/lib
#      "[HKEY_LOCAL_MACHINE\\SOFTWARE\\DRACO\\Current;BinPath]/lib"
    DOC "Path to the DRACO <component> library."
    )

  mark_as_advanced( 
     DRACO_${component}_INCLUDE_DIR 
     DRACO_${component}_LIBRARY
     )

  if(DRACO_${component}_INCLUDE_DIR AND DRACO_${component}_LIBRARY)
    set(DRACO_${component}_FOUND TRUE PARENT_SCOPE)

    list(APPEND DRACO_INCLUDE_DIRS
      ${DRACO_${component}_INCLUDE_DIR}
      )
    list(REMOVE_DUPLICATES DRACO_INCLUDE_DIRS)
    set(DRACO_INCLUDE_DIRS ${DRACO_INCLUDE_DIRS} PARENT_SCOPE)

    list(APPEND DRACO_LIBRARIES
      ${DRACO_${component}_LIBRARY}
      )
    set(DRACO_LIBRARIES ${DRACO_LIBRARIES} PARENT_SCOPE)
    if( NOT DRACO_FIND_QUIETLY )
       message( STATUS "Looking for ${component}...\t${DRACO_${component}_LIBRARY}")
    endif()

 else()

    if( NOT DRACO_FIND_QUIETLY )
       message( STATUS "Looking for ${component}...\t NOT FOUND!\n")
    endif()

  endif()
endfunction(FIND_DRACO_API)

# function(FIND_DRACO_EXE component)
#   set(_DRACO_EXECUTABLE
#     ${DRACO_EXECUTABLE_DIR}/${component}${CMAKE_EXECUTABLE_SUFFIX})
#   if(EXISTS ${_DRACO_EXECUTABLE})
#     set(DRACO_${component}_EXECUTABLE
#       ${_DRACO_EXECUTABLE}
#        PARENT_SCOPE
#        )
#     set(DRACO_${component}_FOUND TRUE PARENT_SCOPE)
#   else(EXISTS ${_DRACO_EXECUTABLE})
#     set(DRACO_${component}_FOUND FALSE PARENT_SCOPE)
#   endif(EXISTS ${_DRACO_EXECUTABLE})
# endfunction(FIND_DRACO_EXE)

#---------------------------------------------------------------------
# Start Actual Work
#---------------------------------------------------------------------

if( NOT DRACO_FIND_QUIETLY )
   message("Looking for DRACO (at $DRACO_DIR and $CMAKE_INSTALL_PREFIX)...")
endif()

# Check for hints found in the environment.
if( "${DRACO_DIR}none" STREQUAL "none" AND EXISTS $ENV{DRACO_DIR} )
   set( DRACO_DIR $ENV{DRACO_DIR} )
endif()
if( "${DRACO_LIBRARY_DIR}none" STREQUAL "none" AND EXISTS $ENV{DRACO_LIB_DIR} )
   set( DRACO_LIBRARY_DIR $ENV{DRACO_LIB_DIR} )
endif()
if( "${DRACO_INCLUDE_DIR}none" STREQUAL "none" AND EXISTS $ENV{DRACO_INCLUDE_DIR} )
   set( DRACO_INCLUDE_DIR $ENV{DRACO_INCLUDE_DIR} )
endif()
if( "${DRACO_INCLUDE_DIR}none" STREQUAL "none" AND EXISTS $ENV{DRACO_INC_DIR} )
   set( DRACO_INCLUDE_DIR $ENV{DRACO_INC_DIR} )
endif()

# Try to find a DRACO installation binary path.
# find_path(DRACO_EXECUTABLE_DIR
#   NAMES mogrify${CMAKE_EXECUTABLE_SUFFIX}
#   PATHS
#     "[HKEY_LOCAL_MACHINE\\SOFTWARE\\DRACO\\Current;BinPath]"
#   DOC "Path to the DRACO binary directory."
#   NO_DEFAULT_PATH
#   )
# find_path(DRACO_EXECUTABLE_DIR
#   NAMES mogrify${CMAKE_EXECUTABLE_SUFFIX}
#   )

# Find each component. Search for all tools in same dir
# <DRACO_EXECUTABLE_DIR>; otherwise they should be found
# independently and not in a cohesive module such as this one.
set(DRACO_FOUND TRUE)
foreach(component ${DRACO_FIND_COMPONENTS} ) # ds++, c4, etc.
   FIND_DRACO_API(${component} ${componenet}/Release.hh )
  # if(component STREQUAL "ds++")
  #   FIND_DRACO_API(ds++ ds++/Release.hh ds++ CORE_RL_ds++_ )
  # elseif(component STREQUAL "MagickWand")
  #   FIND_DRACO_API(MagickWand wand/MagickWand.h
  #     Wand MagickWand CORE_RL_wand_
  #     )
  # elseif(component STREQUAL "MagickCore")
  #   FIND_DRACO_API(MagickCore magick/MagickCore.h
  #     Magick MagickCore CORE_RL_magick_
  #     )
  # else(component STREQUAL "Magick++")
  #   if(DRACO_EXECUTABLE_DIR)
  #     FIND_DRACO_EXE(${component})
  #   endif(DRACO_EXECUTABLE_DIR)
  # endif()
  
  if(NOT DRACO_${component}_FOUND)
    list(FIND DRACO_FIND_COMPONENTS ${component} is_requested)
    if(is_requested GREATER -1)
      set(DRACO_FOUND FALSE)
    endif(is_requested GREATER -1)
  endif(NOT DRACO_${component}_FOUND)
endforeach(component)

# not sure about these
set(DRACO_INCLUDE_DIRS ${DRACO_INCLUDE_DIRS})
set(DRACO_LIBRARIES ${DRACO_LIBRARIES})

get_filename_component( tmp ${DRACO_ds++_LIBRARY} PATH )
if( EXISTS ${tmp} )
   get_filename_component( DRACO_DIR ${tmp} PATH CACHE )
else()
   set( DRACO_DIR ${DRACO_DIR} CACHE FILEPATH 
      "Location of installed DRACO files." )
   message( FATAL_ERROR 
"Could not find ds++ library.  Unable to continue without valid "
"draco.  Please set DRACO_DIR to the Draco installation director !")
endif()

set( DRACO_DIR ${DRACO_DIR} CACHE FILEPATH 
   "Location of installed DRACO files." )
set( DRACO_LIBRARY_DIR ${DRACO_DIR}/lib CACHE FILEPATH 
   "Location of installed DRACO libraries." )
set( DRACO_INCLUDE_DIR ${DRACO_DIR}/include CACHE FILEPATH 
   "Location of installed DRACO libraries." )

#---------------------------------------------------------------------
# Standard Package Output
#---------------------------------------------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args( DRACO DEFAULT_MSG DRACO_FOUND )
# Maintain consistency with all other variables.
set(DRACO_FOUND ${DRACO_FOUND})

#---------------------------------------------------------------------
# DEPRECATED: Setting variables for backward compatibility.
#---------------------------------------------------------------------
# set(DRACO_BINARY_PATH          ${DRACO_EXECUTABLE_DIR}
#     CACHE PATH "Path to the DRACO binary directory.")
# set(DRACO_CONVERT_EXECUTABLE   ${DRACO_convert_EXECUTABLE}
#     CACHE FILEPATH "Path to DRACO's convert executable.")
# set(DRACO_MOGRIFY_EXECUTABLE   ${DRACO_mogrify_EXECUTABLE}
#     CACHE FILEPATH "Path to DRACO's mogrify executable.")
# set(DRACO_IMPORT_EXECUTABLE    ${DRACO_import_EXECUTABLE}
#     CACHE FILEPATH "Path to DRACO's import executable.")
# set(DRACO_MONTAGE_EXECUTABLE   ${DRACO_montage_EXECUTABLE}
#     CACHE FILEPATH "Path to DRACO's montage executable.")
# set(DRACO_COMPOSITE_EXECUTABLE ${DRACO_composite_EXECUTABLE}
#     CACHE FILEPATH "Path to DRACO's composite executable.")
#mark_as_advanced(
#  DRACO_BINARY_PATH
#  DRACO_CONVERT_EXECUTABLE
#  DRACO_MOGRIFY_EXECUTABLE
#  DRACO_IMPORT_EXECUTABLE
#  DRACO_MONTAGE_EXECUTABLE
#  DRACO_COMPOSITE_EXECUTABLE
#  )
