# - Find clapack
# Find the native clapack includes and library
#
#  clapack_INCLUDE_DIRS   - where to find clapack.h, etc.
#  clapack_LIBRARIES      - List of libraries when using clapack.
#  clapack_FOUND          - True if clapack found.

FIND_PATH( LAPACK_INCLUDE_DIR clapack.h
    #"[HKEY_LOCAL_MACHINE\\SOFTWARE\\GnuWin32\\clapack;InstallPath]/include"
    ${LAPACK_INC_DIR}
    $ENV{LAPACK_INC_DIR}
    ${VENDOR_DIR}/atlas/include
)

if( WIN32 )
   set(LAPACK_LIBRARY_NAME lapack.lib)
else()
   set(LAPACK_LIBRARY_NAME lapack)
endif()
FIND_LIBRARY(LAPACK_LIBRARY
    NAMES ${LAPACK_LIBRARY_NAME}
    PATHS
        ${LAPACK_LIB_DIR}
        $ENV{LAPACK_LIB_DIR}
        ${VENDOR_DIR}/atlas/lib
)
if( WIN32 )
   set(BLAS_LIBRARY_NAME blas.lib)
else()
   set(BLAS_LIBRARY_NAME cblas)
endif()

FIND_LIBRARY(LAPACK_BLAS_LIBRARY
    NAMES ${BLAS_LIBRARY_NAME}
    PATHS
        ${LAPACK_LIB_DIR}
        $ENV{LAPACK_LIB_DIR}
        ${VENDOR_DIR}/atlas/lib
)
FIND_LIBRARY(LAPACK_F2C_LIBRARY
    NAMES libf2c.lib
    PATHS
        ${LAPACK_LIB_DIR}
        $ENV{LAPACK_LIB_DIR}
        ${VENDOR_DIR}/atlas/lib
)
FIND_LIBRARY(LAPACK_F77BLAS_LIBRARY
    NAMES f77blas
    PATHS
        ${LAPACK_LIB_DIR}
        $ENV{LAPACK_LIB_DIR}
        ${VENDOR_DIR}/atlas/lib
)
FIND_LIBRARY(LAPACK_ATLAS_LIBRARY
    NAMES atlas
    PATHS
        ${LAPACK_LIB_DIR}
        $ENV{LAPACK_LIB_DIR}
        ${VENDOR_DIR}/atlas/lib
)
MARK_AS_ADVANCED( 
   LAPACK_LIBRARY
   LAPACK_INCLUDE_DIR
   LAPACK_BLAS_LIBRARY
   LAPACK_F2C_LIBRARY 
   LAPACK_F77BLAS_LIBRARY 
   LAPACK_ATLAS_LIBRARY
   )

# handle the QUIETLY and REQUIRED arguments and set clapack_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(
   LAPACK 
   DEFAULT_MSG 
   LAPACK_INCLUDE_DIR
   LAPACK_LIBRARY)

if (LAPACK_FOUND)
    set( LAPACK_FOUND ${LAPACK_FOUND} CACHE BOOL 
       "Did we find the LAPACK libraries?" )
    set(LAPACK_INCLUDE_DIRS ${LAPACK_INCLUDE_DIR})
    set(LAPACK_LIBRARIES    ${LAPACK_LIBRARY} ${LAPACK_BLAS_LIBRARY} )
    if (LAPACK_F2C_LIBRARY)
       list(APPEND LAPACK_LIBRARIES ${LAPACK_F2C_LIBRARY})
    endif()
    if (LAPACK_F77BLAS_LIBRARY)
       list(APPEND LAPACK_LIBRARIES ${LAPACK_F77BLAS_LIBRARY})
    endif()
    if (LAPACK_ATLAS_LIBRARY)
       list(APPEND LAPACK_LIBRARIES ${LAPACK_ATLAS_LIBRARY})
    endif()
endif()

if( VERBOSE )
message("
LAPACK_FOUND =        ${LAPACK_FOUND}
LAPACK_INCLUDE_DIRS = ${LAPACK_INCLUDE_DIR}
LAPACK_LIBRARIES    = ${LAPACK_LIBRARIES}
")
endif()
