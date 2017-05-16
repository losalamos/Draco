#-----------------------------*-cmake-*----------------------------------------#
# file   config/FindEOSPAC.cmake
# date   2017 February 28
# brief  Instructions for discovering the EOSPAC vendor libraries.
# note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# - Find EOSPAC
# Find the native EOSPAC includes and library
#
#  EOSPAC_INCLUDE_DIRS   - where to find eos_Interface.h, etc.
#  EOSPAC_LIBRARIES      - List of libraries when using EOSPAC.
#  EOSPAC_FOUND          - True if EOSPAC found.

if (APPLE)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;.dylib")
endif()

find_path( EOSPAC_INCLUDE_DIR 
    NAMES
       eos_Interface.h
    PATHS
       ${EOSPAC_INC_DIR}
       $ENV{EOSPAC_INC_DIR}
       $ENV{VENDOR_DIR}/include
       ${VENDOR_DIR}/include
    NO_DEFAULT_PATH
)

if( WIN32 )
   if( EOSPAC_STATIC )
      set( EOSPAC_LIBRARY_NAME eospac.lib)
   else()
      set( EOSPAC_LIBRARY_NAME libeospac_dll.lib)
   endif()
else()
   set( EOSPAC_LIBRARY_NAME eospac6)
endif()

find_library(EOSPAC_LIBRARY
    NAMES ${EOSPAC_LIBRARY_NAME}
    PATHS
        ${EOSPAC_LIB_DIR}
        $ENV{EOSPAC_LIB_DIR}
        $ENV{VENDOR_DIR}/lib
        ${VENDOR_DIR}/lib
    NO_DEFAULT_PATH
)

# If above fails, look in default locations
if( NOT EOSPAC_LIBRARY )
   find_path( EOSPAC_INCLUDE_DIR    NAMES eos_Interface.h )
   find_library(EOSPAC_LIBRARY      NAMES ${EOSPAC_LIBRARY_NAME} )
endif()
mark_as_advanced( EOSPAC_LIBRARY EOSPAC_INCLUDE_DIR )

# handle the QUIETLY and REQUIRED arguments and set EOSPAC_FOUND to TRUE if 
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(EOSPAC DEFAULT_MSG EOSPAC_INCLUDE_DIR EOSPAC_LIBRARY)

if (EOSPAC_FOUND)
   set(EOSPAC_FOUND ${EOSPAC_FOUND} CACHE BOOL "Did we find the EOSPAC libraries?")
   set(EOSPAC_INCLUDE_DIRS ${EOSPAC_INCLUDE_DIR})
   set(EOSPAC_LIBRARIES    ${EOSPAC_LIBRARY} CACHE
      FILEPATH "EOSPAC libraries for linking."  )
   
   string( REPLACE "_dll.lib" ".dll" EOSPAC_DLL ${EOSPAC_LIBRARY} )
   mark_as_advanced( EOSPAC_DLL )
   if( EXISTS ${EOSPAC_DLL} )
      set(EOSPAC_DLL_LIBRARIES "${EOSPAC_DLL}" CACHE STRING 
         "list of eospac dll files.")
   mark_as_advanced( EOSPAC_DLL_LIBRARIES EOSPAC_LIBRARIES )
   else()
      set( EOSPAC_DLL "NOTFOUND")
   endif()
endif()

if( VERBOSE )
message("
EOSPAC_FOUND =        ${EOSPAC_FOUND}
EOSPAC_INCLUDE_DIRS = ${EOSPAC_INCLUDE_DIR}
EOSPAC_LIBRARIES    = ${EOSPAC_LIBRARY}
")
endif()
