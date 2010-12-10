# File: dracoVersion.cmake

# Returns:
# DRACO_VERSION_MAJOR
# DRACO_VERSION_MINOR
# DRACO_VERSION_MINOR
# DRACO_VERSION_PATCH

# ${DRACO_DATE_STAMP_YEAR}
# ${DRACO_DATE_STAMP_MONTH}
# ${DRACO_DATE_STAMP_DAY}

set(DRACO_VERSION_MAJOR 6)
set(DRACO_VERSION_MINOR 0)
set(DRACO_VERSION_PATCH 0)

if( UNIX )
    execute_process( COMMAND date +%m/%d/%Y
       OUTPUT_VARIABLE configureDate )
endif()
	
if( WIN32 )
  execute_process( COMMAND "cmd" "/c" "date" "/t"
	OUTPUT_VARIABLE configureDate )
endif()

  string( REGEX REPLACE 
     ".*([0-9][0-9])/([0-9][0-9])/(20[0-9][0-9]).*"
     "\\3" DRACO_DATE_STAMP_YEAR "${configureDate}" )
  string( REGEX REPLACE 
     ".*([0-9][0-9])/([0-9][0-9])/(20[0-9][0-9]).*"
     "\\1" DRACO_DATE_STAMP_MONTH "${configureDate}" )
  string( REGEX REPLACE 
     ".*([0-9][0-9])/([0-9][0-9])/(20[0-9][0-9]).*"
     "\\2" DRACO_DATE_STAMP_DAY "${configureDate}" )

if( "${DRACO_VERSION_MINOR}" MATCHES "[1-9]?[1-9]$")
   set( DRACO_VERSION_PATCH
      "${DRACO_DATE_STAMP_YEAR}${DRACO_DATE_STAMP_MONTH}${DRACO_DATE_STAMP_DAY}"
      )
endif()

set( DRACO_VERSION "${DRACO_VERSION_MAJOR}.${DRACO_VERSION_MINOR}"
   CACHE STRING "Draco version information" FORCE)
set( DRACO_VERSION_FULL  "${DRACO_VERSION}.${DRACO_VERSION_PATCH}"
   CACHE STRING "Draco version information" FORCE)

message( STATUS "This is Draco version ${DRACO_VERSION_FULL}.")

