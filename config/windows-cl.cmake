#-----------------------------*-cmake-*----------------------------------------#
# file   windows-cl.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 June 5
# brief  Establish flags for Windows - MSVC
# note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

#
# Sanity Checks
#
if( NOT ${CMAKE_GENERATOR} MATCHES "Visual Studio" AND
    NOT ${CMAKE_GENERATOR} MATCHES  "NMake Makefiles" )
  message( FATAL_ERROR
  "config/windows-cl.cmake must be taught to build for this compiler "
  "(CMAKE_GENERATOR = ${CMAKE_GENERATOR}). Yell at kt for help on this error." )
endif()

#
# Compiler flag checks
#
include(platform_checks)
query_openmp_availability()

# This is required to provide compatibility between MSVC and MinGW generated libraries.
if( DRACO_SHARED_LIBS )
  set( CMAKE_GNUtoMS ON CACHE BOOL "Compatibility flag for MinGW/MSVC." FORCE)
endif()

# Extra setup (ds++/config.h) for MSVC
# 1. Allow M_PI to be found via <cmath>
set( _USE_MATH_DEFINES 1 )

# Automatically export all symbols. This will add the
# WINDOWS_EXPORT_ALL_SYMBOLS=TRUE target property to all libraries.
# Ref: https://blog.kitware.com/create-dlls-on-windows-without-declspec-using-new-cmake-export-all-feature/
# set( CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE CACHE BOOL "export all library
# symbols." FORCE )

set( MD_or_MT_debug "${MD_or_MT}d" )
if( "${DEBUG_RUNTIME_EXT}" STREQUAL "d" )
  set( MD_or_MT_debug "${MD_or_MT}${DEBUG_RUNTIME_EXT} /RTC1" )
endif()

set( numproc $ENV{NUMBER_OF_PROCESSORS} )
if( "${numproc}notfound" STREQUAL "notfound" )
  set( numproc 1 )
endif()

if( NOT CXX_FLAGS_INITIALIZED )
  set( CXX_FLAGS_INITIALIZED "yes" CACHE INTERNAL "using draco settings." )

  # Notes on options:
  # - /wd 4251 disable warning #4251: 'identifier' : class 'type' needs to have
  #   dll-interface to be used by clients of class 'type2'
  # - /arch:AVX2 At least for VS2013 with NMake Makefiles, CMake doesn't
  #   populate the correct project property. Also, this option causes 'illegal
  #   instruction' for rng on KT's desktop (2017-02-14).
  set( CMAKE_C_FLAGS "/W2 /Gy /fp:precise /DWIN32 /D_WINDOWS /MP${numproc} /wd4251" )
  set( CMAKE_C_FLAGS_DEBUG "/${MD_or_MT_debug} /Od /Zi /DDEBUG /D_DEBUG" )
  set( CMAKE_C_FLAGS_RELEASE "/${MD_or_MT} /O2 /DNDEBUG" )
  set( CMAKE_C_FLAGS_MINSIZEREL "/${MD_or_MT} /O1 /DNDEBUG" )
  set( CMAKE_C_FLAGS_RELWITHDEBINFO "/${MD_or_MT} /O2 /Zi /DDEBUG" )

  # Suppress some MSVC warnings about "unsafe" pointer use.
  if(MSVC_VERSION GREATER 1399)
    string( APPEND CMAKE_C_FLAGS
      " /D_CRT_SECURE_NO_DEPRECATE /D_SCL_SECURE_NO_DEPRECATE /D_SECURE_SCL=0" )
    #string( APPEND CMAKE_C_FLAGS_DEBUG
    #        "${CMAKE_C_FLAGS_DEBUG} /D_HAS_ITERATOR_DEBUGGING=0" )
  endif()

  # If building static libraries, inlcude debugging information in the library.
  if( ${DRACO_LIBRARY_TYPE} MATCHES "STATIC" )
    string( APPEND CMAKE_C_FLAGS_DEBUG " /Z7"   )
  endif()

  set( CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} /EHa" )
  set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}" )
  set( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}" )
  set( CMAKE_CXX_FLAGS_MINSIZEREL "/${MD_or_MT} /O1 /DNDEBUG" )
  set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "/${MD_or_MT} /O2 /Zi /DDEBUG" )

endif()

##---------------------------------------------------------------------------##
# Ensure cache values always match current selection
##---------------------------------------------------------------------------##
set( CMAKE_C_FLAGS                "${CMAKE_C_FLAGS}"                CACHE STRING "compiler flags" FORCE )
set( CMAKE_C_FLAGS_DEBUG          "${CMAKE_C_FLAGS_DEBUG}"          CACHE STRING "compiler flags" FORCE )
set( CMAKE_C_FLAGS_RELEASE        "${CMAKE_C_FLAGS_RELEASE}"        CACHE STRING "compiler flags" FORCE )
set( CMAKE_C_FLAGS_MINSIZEREL     "${CMAKE_C_FLAGS_MINSIZEREL}"     CACHE STRING "compiler flags" FORCE )
set( CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}" CACHE STRING "compiler flags" FORCE )

set( CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS}"                CACHE STRING "compiler flags" FORCE )
set( CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG}"          CACHE STRING "compiler flags" FORCE )
set( CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE}"        CACHE STRING "compiler flags" FORCE )
set( CMAKE_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_MINSIZEREL}"     CACHE STRING "compiler flags" FORCE )
set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}" CACHE STRING "compiler flags" FORCE )

#
# Toggle compiler flags for optional features
#
if( NOT "${OpenMP_C_FLAGS}x" STREQUAL "x" )
  toggle_compiler_flag( OPENMP_FOUND ${OpenMP_C_FLAGS} "C;CXX" "" )
endif()

#
# Extra runtime libraries...
#

find_library( Lib_win_winsock NAMES wsock32;winsock32;ws2_32 )
if( EXISTS "${Lib_win_winsock}" AND CMAKE_CL_64 )
  string(REPLACE "um/x86" "um/x64" Lib_win64_winsock "${Lib_win_winsock}" )
  if( EXISTS "${Lib_win64_winsock}" )
    set( Lib_win_winsock "${Lib_win64_winsock}")
  endif()
endif()

if( ${Lib_win_winsock} MATCHES "NOTFOUND" )
  message( FATAL_ERROR "Could not find library winsock32 or ws2_32!" )
endif()

#------------------------------------------------------------------------------#
# End windows-cl.cmake
#------------------------------------------------------------------------------#
