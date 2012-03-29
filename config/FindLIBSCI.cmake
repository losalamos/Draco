# - Find LIBSCI
# Find the Cray LIBSCI includes and library
#
#  LIBSCI_INCLUDE_DIRS   - where to find LIBSCI.h, etc.
#  LIBSCI_LIBRARIES      - List of libraries when using LIBSCI.
#  LIBSCI_FOUND          - True if LIBSCI found.

if (APPLE)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;.dylib")
endif()

find_path( LIBSCI_INCLUDE_DIR 
    NAMES
       case.h
    PATHS
       ${LIBSCI_INC_DIR}
       $ENV{LIBSCI_INC_DIR}
       $ENV{VENDOR_DIR}/include
       ${VENDOR_DIR}/include
       $ENV{LIBSCI_BASE_DIR}/pgi/include  # Cielito/Cielo
       $ENV{LIBSCI_BASE_DIR}/pgi/119/interlagos/include  # Cielito/Cielo
    NO_DEFAULT_PATH
)

if( WIN32 )
   set( LIBSCI_LIBRARY_NAME sci_dll.lib )
else()
   set( LIBSCI_LIBRARY_NAME sci sci_pgi)
endif()

find_library(LIBSCI_LIBRARY
    NAMES ${LIBSCI_LIBRARY_NAME}
    PATHS
        ${LIBSCI_LIB_DIR}
        $ENV{LIBSCI_LIB_DIR}
        $ENV{VENDOR_DIR}/lib
        ${VENDOR_DIR}/lib
        $ENV{LIBSCI_BASE_DIR}/pgi/lib # Cielito/Cielo
        $ENV{LIBSCI_BASE_DIR}/pgi/119/interlagos/lib # Cielito/Cielo
    NO_DEFAULT_PATH # avoid picking up /usr/lib/liblibsci.so
)

mark_as_advanced( LIBSCI_LIBRARY LIBSCI_INCLUDE_DIR )

# handle the QUIETLY and REQUIRED arguments and set LIBSCI_FOUND to TRUE if 
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LIBSCI DEFAULT_MSG LIBSCI_INCLUDE_DIR LIBSCI_LIBRARY)

if (LIBSCI_FOUND)
   set(LIBSCI_FOUND ${LIBSCI_FOUND} CACHE BOOL "Did we find the LIBSCI libraries?")
   set(LIBSCI_INCLUDE_DIRS ${LIBSCI_INCLUDE_DIR})
   set(LIBSCI_LIBRARIES    ${LIBSCI_LIBRARY} ${LIBSCI_BLAS_LIBRARY} CACHE
      FILEPATH "LIBSCI libraries for linking."  )
   
   string( REPLACE "_dll.lib" ".dll" LIBSCI_DLL ${LIBSCI_LIBRARY} )
   mark_as_advanced( LIBSCI_DLL LIBSCI_BLAS_DLL )
   if( EXISTS ${LIBSCI_DLL} )
      set(LIBSCI_DLL_LIBRARIES "${LIBSCI_DLL};${LIBSCI_BLAS_DLL}" CACHE STRING 
         "list of libsci dll files.")
   mark_as_advanced( LIBSCI_DLL_LIBRARIES LIBSCI_LIBRARIES )
   else()
      set( LIBSCI_DLL "NOTFOUND")
      set( LIBSCI_BLAS_DLL "NOTFOUND" )
   endif()
endif()

# libsci provides the following tools:

set( BLAS_FOUND       ${LIBSCI_FOUND} )
set( BLAS_LIBRARIES   ${LIBSCI_LIBRARY} )
set( LAPACK_FOUND     ${LIBSCI_FOUND} )
set( LAPACK_LIBRARIES ${LIBSCI_LIBRARY} )

if( VERBOSE )
message("
LIBSCI_FOUND =        ${LIBSCI_FOUND}
LIBSCI_INCLUDE_DIRS = ${LIBSCI_INCLUDE_DIR}
LIBSCI_LIBRARIES    = ${LIBSCI_LIBRARY}
BLAS_LIBRARIES      = ${BLAS_LIBRARIES}
LAPACK_LIBRARIES    = ${LAPACK_LIBRARIES}
")
endif()
