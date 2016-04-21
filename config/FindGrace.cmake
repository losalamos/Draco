#-----------------------------*-cmake-*----------------------------------------#
# file   config/FindGrace.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   2011 September 29
# brief  Instructions for discovering the Grace vendor libraries.
# note   Copyright (C) 2016 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# - Find Grace
# Find the native Grace includes and library
#
#  Grace_INCLUDE_DIRS   - where to find Grace.h, etc.
#  Grace_LIBRARIES      - List of libraries when using Grace.
#  GRACE_FOUND          - True if Grace found.

if (APPLE)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;.dylib")
endif()

find_path( Grace_INCLUDE_DIR
    NAMES
       grace_np.h
    PATHS
       ${GRACE_INC_DIR}
       $ENV{GRACE_INC_DIR}
       $ENV{VENDOR_DIR}/include
       ${VENDOR_DIR}/include
    NO_DEFAULT_PATH
)

if( WIN32 )
   set( Grace_LIBRARY_NAME libgrace_dll.lib)
else()
   set( Grace_LIBRARY_NAME grace_np)
endif()

find_library(Grace_LIBRARY
    NAMES ${Grace_LIBRARY_NAME}
    PATHS
        ${GRACE_LIB_DIR}
        $ENV{GRACE_LIB_DIR}
        $ENV{VENDOR_DIR}/lib
        ${VENDOR_DIR}/lib
    NO_DEFAULT_PATH
)

mark_as_advanced( Grace_LIBRARY Grace_INCLUDE_DIR )

# handle the QUIETLY and REQUIRED arguments and set GRACE_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Grace DEFAULT_MSG Grace_INCLUDE_DIR Grace_LIBRARY)

if( GRACE_FOUND )
   set(GRACE_FOUND ${GRACE_FOUND} CACHE BOOL "Did we find the Grace libraries?")
   set(Grace_INCLUDE_DIRS ${Grace_INCLUDE_DIR})
   set(Grace_LIBRARIES    ${Grace_LIBRARY} ${Grace_BLAS_LIBRARY} CACHE
      FILEPATH "Grace libraries for linking."  )
   mark_as_advanced( Grace_LIBRARIES )

endif()

if( VERBOSE )
message("
GRACE_FOUND =        ${GRACE_FOUND}
Grace_INCLUDE_DIRS = ${Grace_INCLUDE_DIR}
Grace_LIBRARIES    = ${Grace_LIBRARY}
")
endif()
