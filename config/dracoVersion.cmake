# File: dracoVersion.cmake

# Returns:
# ${PROJNAME}_VERSION_MAJOR
# ${PROJNAME}_VERSION_MINOR
# ${PROJNAME}_VERSION_MINOR
# ${PROJNAME}_VERSION_PATCH

# ${${PROJNAME}_DATE_STAMP_YEAR}
# ${${PROJNAME}_DATE_STAMP_MONTH}
# ${${PROJNAME}_DATE_STAMP_DAY}

macro( set_ccs2_software_version PROJNAME )

   if( "${${PROJNAME}_VERSION_MAJOR}notset" STREQUAL "notset" )
      message( WARNING "${PROJNAME}_VERSION_MAJOR should already be set!" )
      set(${PROJNAME}_VERSION_MAJOR 0)
   endif()
   if( "${${PROJNAME}_VERSION_MINOR}notset" STREQUAL "notset" )
      message( WARNING "${PROJNAME}_VERSION_MINOR should already be set!" )
      set(${PROJNAME}_VERSION_MINOR 0)
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
      "\\3" ${PROJNAME}_DATE_STAMP_YEAR "${configureDate}" )
   string( REGEX REPLACE 
      ".*([0-9][0-9])/([0-9][0-9])/(20[0-9][0-9]).*"
      "\\1" ${PROJNAME}_DATE_STAMP_MONTH "${configureDate}" )
   string( REGEX REPLACE 
      ".*([0-9][0-9])/([0-9][0-9])/(20[0-9][0-9]).*"
      "\\2" ${PROJNAME}_DATE_STAMP_DAY "${configureDate}" )
   
   if( "${${PROJNAME}_VERSION_PATCH}notset" STREQUAL "notset" ) # "[1-9]?[1-9]$")
      set( ${PROJNAME}_VERSION_PATCH
         "${${PROJNAME}_DATE_STAMP_YEAR}${${PROJNAME}_DATE_STAMP_MONTH}${${PROJNAME}_DATE_STAMP_DAY}"
         )
   endif()
   
   set( ${PROJNAME}_VERSION "${${PROJNAME}_VERSION_MAJOR}.${${PROJNAME}_VERSION_MINOR}"
      CACHE STRING "${PROJNAME} version information" FORCE)
   set( ${PROJNAME}_VERSION_FULL  "${${PROJNAME}_VERSION}.${${PROJNAME}_VERSION_PATCH}"
      CACHE STRING "${PROJNAME} version information" FORCE)
   
   message( STATUS "This is ${PROJNAME} version ${${PROJNAME}_VERSION_FULL}.")
   
   
   # Support for CPack
   set( CPACK_PACKAGE_VERSION_MAJOR ${${PROJNAME}_VERSION_MAJOR} )
   set( CPACK_PACKAGE_VERSION_MINOR ${${PROJNAME}_VERSION_MINOR} )
   set( CPACK_PACKAGE_VERSION_PATCH ${${PROJNAME}_VERSION_PATCH} )
   
endmacro()