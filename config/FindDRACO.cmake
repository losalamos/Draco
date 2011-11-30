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
#---------------------------------------------------------------------


#---------------------------------------------------------------------
# PADSTRING( result ${length} ${origstring} )
#---------------------------------------------------------------------
macro(PADSTRING length origstring )

   # copy the orig string to the result.
   set( paddedstring ${origstring} )
   string( LENGTH ${origstring} len_of_origstring )

   # Loop over the remainder of the desired lenth, appending a space
   # for each char.
   set( index ${len_of_origstring} )
   while( ${index} LESS ${length} )
      set( paddedstring "${paddedstring} " )
      math( EXPR index '${index}+1' )
   endwhile()

endmacro()

#---------------------------------------------------------------------
# FIND_DRACO_API(ds++ SP.hh ds++ CORE_RL_ds++_ )
#---------------------------------------------------------------------
function(FIND_DRACO_API component header)
  set(DRACO_${component}_FOUND FALSE PARENT_SCOPE)

  find_path(DRACO_${component}_INCLUDE_DIR
    NAMES ${component}
    PATHS
      ${DRACO_INCLUDE_DIR}
      ${DRACO_DIR}/include
      ${CMAKE_INSTALL_PREFIX}/include
#    PATH_SUFFIXES
#       ${component}
    DOC "Path to the DRACO <component> include dir."
    )
  find_library(DRACO_${component}_LIBRARY
    NAMES rtt_${component}
    PATHS
      ${DRACO_LIBRARY_DIR}
      ${DRACO_DIR}/lib
      ${CMAKE_INSTALL_PREFIX}/lib
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
       set(length 20)
       PADSTRING( ${length} ${component} ) # returns ${paddedstring}
       message( STATUS "Looking for ${paddedstring}...\t${DRACO_${component}_LIBRARY}")
    endif()

 else()

    if( NOT DRACO_FIND_QUIETLY )
       message( STATUS "Looking for ${component}...\t NOT FOUND!\n")
    endif()

  endif()
endfunction(FIND_DRACO_API)

#---------------------------------------------------------------------
# Start Actual Work
#---------------------------------------------------------------------

if( NOT DRACO_FIND_QUIETLY )
   message("
Looking for DRACO...
-- Searching DRACO_DIR            = ${DRACO_DIR} 
             CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}")
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

# Find each component. Search for all tools in same dir
# <DRACO_EXECUTABLE_DIR>; otherwise they should be found
# independently and not in a cohesive module such as this one.
set(DRACO_FOUND TRUE)
foreach(component ${DRACO_FIND_COMPONENTS} ) # ds++, c4, etc.
   FIND_DRACO_API(${component} ${componenet}/Release.hh )
  
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


# Find out if DRACO_C4 is SCALAR or MPI
if( DRACO_FOUND )
   if( EXISTS ${DRACO_DIR}/config/Dracoo-CMakeCache.txt )
      # Read complete Draco CMakeCache.txt file and store in variable
      file( READ ${DRACO_DIR}/config/Dracoo-CMakeCache.txt filedata )
      # Convert string into lines by replacing EOL with semi-colon.
      string( REGEX REPLACE "\n" ";" linedata ${filedata} )
      # Loop over lines looking for DRACO_C4 
      foreach( line ${filedata} )
         if( "${line}" MATCHES "DRACO_C4" )
            message("${line}")
         endif()
      endforeach()
   endif()
endif()
