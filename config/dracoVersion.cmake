# File: dracoVersion.cmake

# Returns:
# DRACO_VERSION_MAJOR
# DRACO_VERSION_MINOR
# DRACO_VERSION_MINOR
# DRACO_VERSION_PATCH

# ${DRACO_DATE_STAMP_YEAR}
# ${DRACO_DATE_STAMP_MONTH}
# ${DRACO_DATE_STAMP_DAY}

if( "${DRACO_VERSION_MAJOR}notset" STREQUAL "notset" )
   message( WARNING "DRACO_VERSION_MAJOR should already be set!" )
   set(DRACO_VERSION_MAJOR 0)
endif()
if( "${DRACO_VERSION_MINOR}notset" STREQUAL "notset" )
   message( WARNING "DRACO_VERSION_MINOR should already be set!" )
   set(DRACO_VERSION_MINOR 0)
endif()

# Use the configure date as the patch number for development builds
# (non release builds)
if( UNIX )
   execute_process( 
      COMMAND         date +%m/%d/%Y
      OUTPUT_VARIABLE configureDate )
elseif( WIN32 )
   execute_process( 
      COMMAND         "cmd" "/c" "date" "/t"
      OUTPUT_VARIABLE configureDate )
endif()

# Format the configureDate
string( REGEX REPLACE 
   ".*([0-9][0-9])/([0-9][0-9])/(20[0-9][0-9]).*"
   "\\3" DRACO_DATE_STAMP_YEAR "${configureDate}" )
string( REGEX REPLACE 
   ".*([0-9][0-9])/([0-9][0-9])/(20[0-9][0-9]).*"
   "\\1" DRACO_DATE_STAMP_MONTH "${configureDate}" )
string( REGEX REPLACE 
   ".*([0-9][0-9])/([0-9][0-9])/(20[0-9][0-9]).*"
   "\\2" DRACO_DATE_STAMP_DAY "${configureDate}" )

if( "${DRACO_VERSION_PATCH}notset" STREQUAL "notset" ) # "[1-9]?[1-9]$")
   set( DRACO_VERSION_PATCH
      "${DRACO_DATE_STAMP_YEAR}${DRACO_DATE_STAMP_MONTH}${DRACO_DATE_STAMP_DAY}"
      )
endif()

set( DRACO_VERSION "${DRACO_VERSION_MAJOR}.${DRACO_VERSION_MINOR}"
   CACHE STRING "Draco version information" FORCE)
set( DRACO_VERSION_FULL  "${DRACO_VERSION}.${DRACO_VERSION_PATCH}"
   CACHE STRING "Draco version information" FORCE)

message( STATUS "This is Draco version ${DRACO_VERSION_FULL}.")


# Support for CPack
set( CPACK_PACKAGE_VERSION_MAJOR ${DRACO_VERSION_MAJOR} )
set( CPACK_PACKAGE_VERSION_MINOR ${DRACO_VERSION_MINOR} )
set( CPACK_PACKAGE_VERSION_PATCH ${DRACO_VERSION_PATCH} )
