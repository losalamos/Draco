#-----------------------------*-cmake-*----------------------------------------#
# file   config/dracoVersion.cmake
# author Kelly G. Thompson, kgt@lanl.gov
# date   2010 Dec 1
# brief  Ensure version is set and use config date as ver. patch value.
# note   Copyright (C) 2016-2019 Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

macro( set_ccs2_software_version PROJNAME )

  if( NOT DEFINED ${PROJNAME}_VERSION_MAJOR )
    message( WARNING "${PROJNAME}_VERSION_MAJOR should already be set!" )
    set(${PROJNAME}_VERSION_MAJOR 0)
  endif()
  if( NOT DEFINED ${PROJNAME}_VERSION_MINOR )
    message( WARNING "${PROJNAME}_VERSION_MINOR should already be set!" )
    set(${PROJNAME}_VERSION_MINOR 0)
  endif()

  # Use the configure date as the patch number for development builds (non-
  # release builds)
  if( UNIX )
    execute_process(
      COMMAND         date +%m/%d/%Y
      OUTPUT_VARIABLE configureDate )
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
  elseif( WIN32 )
    execute_process(
      COMMAND         "cmd" "/c" "date" "/t"
      OUTPUT_VARIABLE configureDate
      OUTPUT_STRIP_TRAILING_WHITESPACE )
    # this should produce a string of the form "Sat MM/DD/YYYY".
    # Format the configureDate
    if( ${configureDate} MATCHES ".*(20[0-9][0-9])-([0-9][0-9])-([0-9][0-9]).*" )
      ### YYYY-MM-DD
      string( REGEX REPLACE
	".*(20[0-9][0-9])-([0-9][0-9])-([0-9][0-9]).*"
	"\\1" ${PROJNAME}_DATE_STAMP_YEAR "${configureDate}" )
      string( REGEX REPLACE
	".*(20[0-9][0-9])-([0-9][0-9])-([0-9][0-9]).*"
	"\\2" ${PROJNAME}_DATE_STAMP_MONTH "${configureDate}" )
      string( REGEX REPLACE
	".*(20[0-9][0-9])-([0-9][0-9])-([0-9][0-9]).*"
	"\\3" ${PROJNAME}_DATE_STAMP_DAY "${configureDate}" )
    else()
      # this should produce a string of the form "Sat MM/DD/YYYY".
      string( REGEX REPLACE
	".*([0-9][0-9])/([0-9][0-9])/(20[0-9][0-9]).*"
	"\\3" ${PROJNAME}_DATE_STAMP_YEAR "${configureDate}" )
      string( REGEX REPLACE
	".*([0-9][0-9])/([0-9][0-9])/(20[0-9][0-9]).*"
	"\\1" ${PROJNAME}_DATE_STAMP_MONTH "${configureDate}" )
      string( REGEX REPLACE
	".*([0-9][0-9])/([0-9][0-9])/(20[0-9][0-9]).*"
	"\\2" ${PROJNAME}_DATE_STAMP_DAY "${configureDate}" )
    endif()
    # message("
    # YYYY = ${${PROJNAME}_DATE_STAMP_YEAR}
    # MM   = ${${PROJNAME}_DATE_STAMP_MONTH}
    # DD   = ${${PROJNAME}_DATE_STAMP_DAY}
    # ")
  endif()

  # Query git branch name?
  # git rev-parse --abbrev-ref HEAD

  string(TOUPPER ${PROJNAME} PROJNAME_UPPER)
  set( ${PROJNAME}_BUILD_DATE
    "${${PROJNAME}_DATE_STAMP_YEAR}/${${PROJNAME}_DATE_STAMP_MONTH}/${${PROJNAME}_DATE_STAMP_DAY}" )
  if( DEFINED ${PROJNAME_UPPER}_VERSION_PATCH ) # "[1-9]?[1-9]$")
    set( ${PROJNAME}_VERSION_PATCH ${${PROJNAME_UPPER}_VERSION_PATCH} )
  else()
    set( ${PROJNAME}_VERSION_PATCH
      "${${PROJNAME}_DATE_STAMP_YEAR}${${PROJNAME}_DATE_STAMP_MONTH}${${PROJNAME}_DATE_STAMP_DAY}"
      )
  endif()

  set( ${PROJNAME}_VERSION "${${PROJNAME}_VERSION_MAJOR}.${${PROJNAME}_VERSION_MINOR}"
    CACHE STRING "${PROJNAME} version information" FORCE)
  set( ${PROJNAME}_VERSION_FULL  "${${PROJNAME}_VERSION}.${${PROJNAME}_VERSION_PATCH}"
    CACHE STRING "${PROJNAME} version information" FORCE)
  mark_as_advanced( ${PROJNAME}_VERSION )

  message( "\n======================================================\n"
    "This is ${PROJNAME} version ${${PROJNAME}_VERSION_FULL}.\n"
     "======================================================\n" )

  # Support for CPack
  set( CPACK_PACKAGE_VERSION_MAJOR ${${PROJNAME}_VERSION_MAJOR} )
  set( CPACK_PACKAGE_VERSION_MINOR ${${PROJNAME}_VERSION_MINOR} )
  set( CPACK_PACKAGE_VERSION_PATCH ${${PROJNAME}_VERSION_PATCH} )

endmacro()

#------------------------------------------------------------------------------#
# End config/dracoVersion.cmake
#------------------------------------------------------------------------------#
