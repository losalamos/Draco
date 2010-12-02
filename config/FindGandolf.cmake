# - Find GANDOLF
# Find the native GANDOLF includes and library
#
#  GANDOLF_LIBRARIES      - List of libraries when using GANDOLF.
#  GANDOLF_FOUND          - True if GANDOLF found.

set( GANDOLF_LIBRARY_NAME gandolf)
set( GANDOLF_GFORTRAN_NAME gfortran)

find_library(GANDOLF_LIBRARY
    NAMES ${GANDOLF_LIBRARY_NAME}
    PATHS
        ${GANDOLF_LIB_DIR}
        $ENV{GANDOLF_LIB_DIR}
)
find_library(GANDOLF_GFORTRAN_LIBRARY
    NAMES ${GANDOLF_GFORTRAN_NAME}
    PATHS
        ${GANDOLF_LIB_DIR}
        $ENV{GANDOLF_LIB_DIR}
)

mark_as_advanced( GANDOLF_LIBRARY GANDOLF_GFORTRAN_LIBRARY )

# handle the QUIETLY and REQUIRED arguments and set GANDOLF_FOUND to TRUE if 
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GANDOLF DEFAULT_MSG GANDOLF_LIBRARY)

if (GANDOLF_FOUND)
#    set(GANDOLF_INCLUDE_DIRS ${GANDOLF_INCLUDE_DIR})
    set(GANDOLF_LIBRARIES    ${GANDOLF_LIBRARY} ${GANDOLF_GFORTRAN_LIBRARY} )
    set(GANDOLF_FOUND ${GANDOLF_FOUND} CACHE BOOL 
       "Did we find the gadolf libraries?")
#    string( REPLACE "_dll.lib" ".dll" GANDOLF_DLL ${GANDOLF_LIBRARY} )
#    string( REPLACE "_dll.lib" ".dll" GANDOLF_BLAS_DLL ${GANDOLF_BLAS_LIBRARY} )
#   mark_as_advanced( GANDOLF_DLL GANDOLF_BLAS_DLL )
#    if( EXISTS ${GANDOLF_DLL} )
#      set(GANDOLF_DLL_LIBRARIES "${GANDOLF_DLL};${GANDOLF_BLAS_DLL}" CACHE STRING 
#         "list of gsl dll files.")
#    else()
#      set( GANDOLF_DLL "NOTFOUND")
#      set( GANDOLF_BLAS_DLL "NOTFOUND" )
#    endif()

endif()

if( VERBOSE )
message("
GANDOLF_FOUND =        ${GANDOLF_FOUND}
GANDOLF_LIBRARIES    = ${GANDOLF_LIBRARY}
")
endif()

