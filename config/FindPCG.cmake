# - Find PCG
# Find the native PCG includes and library
#
#  PCG_LIBRARIES      - List of libraries when using PCG.
#  PCG_FOUND          - True if PCG found.

if (APPLE)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;.dylib")
endif()

set( PCG_LIBRARY_NAME pcg)

find_library(PCG_LIBRARY
    NAMES ${PCG_LIBRARY_NAME}
    PATHS
        ${PCG_LIB_DIR}
        $ENV{PCG_LIB_DIR}
        $ENV{VENDOR_DIR}/lib
        ${VENDOR_DIR}/lib
    NO_DEFAULT_PATH # avoid picking up /usr/lib/libgsl.so
)
mark_as_advanced( PCG_LIBRARY PCG_LIBRARY )

# handle the QUIETLY and REQUIRED arguments and set PCG_FOUND to TRUE if 
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PCG DEFAULT_MSG PCG_LIBRARY)

if (PCG_FOUND)
   set(PCG_FOUND ${PCG_FOUND} CACHE BOOL "Did we find the PCG libraries?")
#   set(PCG_INCLUDE_DIRS ${PCG_INCLUDE_DIR})
#   set(PCG_LIBRARIES    ${PCG_LIBRARY} ${PCG_BLAS_LIBRARY} CACHE
#      FILEPATH "PCG libraries for linking."  )
   
#   string( REPLACE "_dll.lib" ".dll" PCG_DLL ${PCG_LIBRARY} )
#   string( REPLACE "_dll.lib" ".dll" PCG_BLAS_DLL ${PCG_BLAS_LIBRARY} )
#   mark_as_advanced( PCG_DLL PCG_BLAS_DLL )
#   if( EXISTS ${PCG_DLL} )
#      set(PCG_DLL_LIBRARIES "${PCG_DLL};${PCG_BLAS_DLL}" CACHE STRING 
#         "list of gsl dll files.")
#   mark_as_advanced( PCG_DLL_LIBRARIES PCG_LIBRARIES )
#   else()
#      set( PCG_DLL "NOTFOUND")
#      set( PCG_BLAS_DLL "NOTFOUND" )
#   endif()
endif()

if( VERBOSE )
   message("
PCG_FOUND   = ${PCG_FOUND}
PCG_LIBRARY = ${PCG_LIBRARY}
")
endif()
