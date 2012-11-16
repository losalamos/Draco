# - Find Random123
# Find the Random123 includes
#
#  RANDOM123_INCLUDE_DIRS   - where to find threefry.h, etc.
#  RANDOM123_FOUND          - True if Random123 found.

find_path( RANDOM123_INCLUDE_DIR
    NAMES
       "Random123/threefry.h"
    PATHS
       ${RANDOM123_INC_DIR}
       $ENV{RANDOM123_INC_DIR}
    NO_DEFAULT_PATH
)

mark_as_advanced( RANDOM123_INCLUDE_DIR )

# handle the QUIETLY and REQUIRED arguments and set RANDOM123_FOUND to
# TRUE if all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Random123 DEFAULT_MSG RANDOM123_INCLUDE_DIR)

if( RANDOM123_FOUND )
   set(RANDOM123_FOUND ${RANDOM123_FOUND} CACHE BOOL "Did we find the Random123 include directory?")
   set(RANDOM123_INCLUDE_DIRS ${RANDOM_INCLUDE_DIR})

   # Enable the definition of UINT64_C in stdint.h
   add_definitions(-D__STDC_CONSTANT_MACROS)
endif()

if( VERBOSE )
message("
RANDOM123_FOUND        = ${RANDOM123_FOUND}
RANDOM123_INCLUDE_DIRS = ${RANDOM123_INCLUDE_DIR}
")
endif()
