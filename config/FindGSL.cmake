# - Find GSL
# Find the native GSL includes and library
#
#  GSL_INCLUDE_DIRS   - where to find GSL.h, etc.
#  GSL_LIBRARIES      - List of libraries when using GSL.
#  GSL_FOUND          - True if GSL found.

find_path( GSL_INCLUDE_DIR 
    NAMES
       gsl/gsl_sf.h
    PATHS
       ${GSL_INC_DIR}
       $ENV{GSL_INC_DIR}
    NO_DEFAULT_PATH
)

if( WIN32 )
   set( GSL_LIBRARY_NAME libgsl_dll.lib)
   set( GSL_BLAS_NAME libgslcblas_dll.lib)
else()
   set( GSL_LIBRARY_NAME gsl)
   set( GSL_BLAS_NAME gslcblas)
endif()

find_library(GSL_LIBRARY
    NAMES ${GSL_LIBRARY_NAME}
    PATHS
        ${GSL_LIB_DIR}
        $ENV{GSL_LIB_DIR}
    NO_DEFAULT_PATH # avoid picking up /usr/lib/libgsl.so
)
find_library(GSL_BLAS_LIBRARY
    NAMES ${GSL_BLAS_NAME}
    PATHS
        ${GSL_LIB_DIR}
        $ENV{GSL_LIB_DIR}
    NO_DEFAULT_PATH
)
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
