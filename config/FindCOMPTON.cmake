#-----------------------------*-cmake-*----------------------------------------#
# file   config/FindCompton.cmake
# author Kendra Keady <keadyk@lanl.gov>
# date   2017 February 28
# brief  Instructions for discovering the Compton vendor libraries.
# note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# - Find Compton
# Find the native Compton includes and library
#
#  COMPTON_INCLUDE_DIRS   - where to find multigroup_lib_builder.hh, etc.
#  COMPTON_LIBRARIES      - List of libraries when using Compton.
#  COMPTON_FOUND          - True if Compton found.

if (APPLE)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;.dylib")
endif()

 find_path( COMPTON_INCLUDE_DIR 
    NAMES
       multigroup_lib_builder.hh
    PATHS
       ${COMPTON_INC_DIR}
       $ENV{COMPTON_INC_DIR}
       $ENV{VENDOR_DIR}/include
       ${VENDOR_DIR}/include
    NO_DEFAULT_PATH
)

# TODO: Win32 logic untested as of 2/28/17.
if( WIN32 )
   if( COMPTON_STATIC )
      set( COMPTON_LIBRARY_NAME libcompton.lib)
   else()
      set( COMPTON_LIBRARY_NAME libcompton_dll.lib)
   endif()
else()
   set( COMPTON_LIBRARY_NAME libLib_compton.a)
endif()

find_library(COMPTON_LIBRARY
    NAMES ${COMPTON_LIBRARY_NAME}
    PATHS
        ${COMPTON_LIB_DIR}
        $ENV{COMPTON_LIB_DIR}
        $ENV{VENDOR_DIR}/lib
        ${VENDOR_DIR}/lib
    NO_DEFAULT_PATH
)

# If above fails, look in default locations
if( NOT COMPTON_LIBRARY )
   find_path( COMPTON_INCLUDE_DIR    NAMES multigroup_lib_builder.hh )
   find_library(COMPTON_LIBRARY      NAMES ${COMPTON_LIBRARY_NAME} )
endif()
mark_as_advanced( COMPTON_LIBRARY COMPTON_INCLUDE_DIR )

# handle the QUIETLY and REQUIRED arguments and set COMPTON_FOUND to TRUE if 
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(COMPTON DEFAULT_MSG COMPTON_INCLUDE_DIR COMPTON_LIBRARY)

if (COMPTON_FOUND)
   set(COMPTON_FOUND ${COMPTON_FOUND} CACHE BOOL "Did we find the Compton libraries?")
   set(COMPTON_INCLUDE_DIRS ${COMPTON_INCLUDE_DIR})
   set(COMPTON_LIBRARIES    ${COMPTON_LIBRARY} CACHE
      FILEPATH "Compton libraries for linking."  )
   
   string( REPLACE "_dll.lib" ".dll" COMPTON_DLL ${COMPTON_LIBRARY} )
   mark_as_advanced( COMPTON_DLL )
   if( EXISTS ${COMPTON_DLL} )
      set(COMPTON_DLL_LIBRARIES "${COMPTON_DLL}" CACHE STRING 
         "list of compton dll files.")
   mark_as_advanced( COMPTON_DLL_LIBRARIES COMPTON_LIBRARIES )
   else()
      set( COMPTON_DLL "NOTFOUND")
   endif()
endif()

if( VERBOSE )
message("
COMPTON_FOUND =        ${COMPTON_FOUND}
COMPTON_INCLUDE_DIRS = ${COMPTON_INCLUDE_DIR}
COMPTON_LIBRARIES    = ${COMPTON_LIBRARY}
")
endif()
