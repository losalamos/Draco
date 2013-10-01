# - Find GSL
# Find the native GSL includes and library
#
#  GSL_INCLUDE_DIRS   - where to find GSL.h, etc.
#  GSL_LIBRARIES      - List of libraries when using GSL.
#  GSL_FOUND          - True if GSL found.

if (APPLE)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;.dylib")
endif()

find_path( GSL_INCLUDE_DIR 
    NAMES
       gsl/gsl_sf.h
    PATHS
       ${GSL_INC_DIR}
       $ENV{GSL_INC_DIR}
       $ENV{VENDOR_DIR}/include
       ${VENDOR_DIR}/include
       $ENV{GSL_INC}/..           # Cielito/Cielo
    NO_DEFAULT_PATH
)

set( GSL_LIBRARY_NAME gsl)
if( WIN32 )
  set( GSL_BLAS_NAME cblas)
  # Both dll and static builds link against the .lib.  However, the dll build
  # must copy the dll to the run directory for tests.
  #if( NOT GSL_STATIC )
    # set(CMAKE_FIND_LIBRARY_PREFIXES "" )
    # set( CMAKE_FIND_LIBRARY_SUFFIXES ".dll" )
  # endif()
else() # Linux
  set( GSL_BLAS_NAME gslcblas)
  if( GSL_STATIC )
    # set(CMAKE_FIND_LIBRARY_PREFIXES lib )
    set( CMAKE_FIND_LIBRARY_SUFFIXES ".a" )
  endif()
endif()

# Use CMake style paths (this conversion is needed on Win32)
file( TO_CMAKE_PATH ${GSL_LIB_DIR} GSL_LIB_DIR )

# For Windows, consider looking for both Release AND Debug versions.  
# If both are found, then GSL_LIBRARY can be set like this:
# set( GSL_LIBRARY
#      Debug ${GSL_LIBRARY_DEBUG}
#      Release ${GSL_LIBRARY_RELEASE} )

find_library(GSL_LIBRARY
    NAMES ${GSL_LIBRARY_NAME}
    PATHS
        ${GSL_LIB_DIR}
        $ENV{GSL_LIB_DIR}
        $ENV{VENDOR_DIR}/lib
        ${VENDOR_DIR}/lib
        $ENV{GSL_DIR}               # Cielito/Cielo
    PATH_SUFFIXES Release Debug
    NO_DEFAULT_PATH                 # avoid picking up /usr/lib/libgsl.so
)

find_library(GSL_BLAS_LIBRARY
    NAMES ${GSL_BLAS_NAME}
    PATHS
        ${GSL_LIB_DIR}
        $ENV{GSL_LIB_DIR}
        $ENV{VENDOR_DIR}/lib
        ${VENDOR_DIR}/lib
        $ENV{GSL_DIR}               # Cielito/Cielo
    PATH_SUFFIXES Release Debug
    NO_DEFAULT_PATH
)
# If above fails, look in default locations
if( NOT GSL_LIBRARY )
   find_path( GSL_INCLUDE_DIR    NAMES gsl/gsl_sf.h )
   find_library(GSL_LIBRARY      NAMES ${GSL_LIBRARY_NAME} )
   find_library(GSL_BLAS_LIBRARY NAMES ${GSL_BLAS_NAME} )
endif()
mark_as_advanced( GSL_LIBRARY GSL_INCLUDE_DIR GSL_LIBRARY GSL_BLAS_LIBRARY )

# handle the QUIETLY and REQUIRED arguments and set GSL_FOUND to TRUE if 
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GSL DEFAULT_MSG GSL_INCLUDE_DIR GSL_LIBRARY)

if (GSL_FOUND)
   set(GSL_FOUND ${GSL_FOUND} CACHE BOOL "Did we find the GSL libraries?")
   set(GSL_INCLUDE_DIRS ${GSL_INCLUDE_DIR})
   set(GSL_LIBRARIES    ${GSL_LIBRARY} ${GSL_BLAS_LIBRARY} CACHE
      FILEPATH "GSL libraries for linking."  )
   
   string( REPLACE "_dll.lib" ".dll" GSL_DLL ${GSL_LIBRARY} )
   string( REPLACE "_dll.lib" ".dll" GSL_BLAS_DLL ${GSL_BLAS_LIBRARY} )
   mark_as_advanced( GSL_DLL GSL_BLAS_DLL )
   if( EXISTS ${GSL_DLL} )
      set(GSL_DLL_LIBRARIES "${GSL_DLL};${GSL_BLAS_DLL}" CACHE STRING 
         "list of gsl dll files.")
   mark_as_advanced( GSL_DLL_LIBRARIES GSL_LIBRARIES )
   else()
      set( GSL_DLL "NOTFOUND")
      set( GSL_BLAS_DLL "NOTFOUND" )
   endif()
endif()

if( VERBOSE )
message("
GSL_FOUND =        ${GSL_FOUND}
GSL_INCLUDE_DIRS = ${GSL_INCLUDE_DIR}
GSL_LIBRARIES    = ${GSL_LIBRARY}
")
endif()
