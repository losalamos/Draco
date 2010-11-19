#-----------------------------*-cmake-*----------------------------------------#
# file   config/locate_draco.cmake
# author 
# date   2010 Nov 2
# brief  Locate draco
# note   Copyright © 2010 LANS, LLC  
#------------------------------------------------------------------------------#
# $Id$ 
#------------------------------------------------------------------------------#

# Try to find Draco.
# 1. Look for installed draco.  User can set DRACO_INC_DIR and
#    DRACO_LIB_DIR which will be searched.  The script will also look
#    at CMAKE_INSTALL_PREFIX. 
# 2. If Draco is not found, then try to download and build it.

message("Looking for Draco...")

# Check for hints found in the environment.
if( "${DRACO_DIR}none" STREQUAL "none" AND 
    EXISTS $ENV{DRACO_DIR} )
   set( DRACO_DIR $ENV{DRACO_DIR} )
endif()
if( "${DRACO_LIB_DIR}none" STREQUAL "none" AND
    EXISTS $ENV{DRACO_LIB_DIR} )
   set( DRACO_LIB_DIR $ENV{DRACO_LIB_DIR} )
endif()
if( "${DRACO_INCLUDE_DIR}none" STREQUAL "none" AND
    EXISTS $ENV{DRACO_INCLUDE_DIR} )
   set( DRACO_INCLUDE_DIR $ENV{DRACO_INCLUDE_DIR} )
endif()

# Find the headers
find_path( DRACO_INCLUDE_DIR ds++/config.h
    ${DRACO_INC_DIR}
    ${DRACO_DIR}/include
    ${CMAKE_INSTALL_PREFIX}/include 
 )
set( components ds++ c4 viz rng mesh_element RTT_Format_Reader
    meshReaders traits cdi cdi_analytic)
foreach( component ${components} )
   string( REGEX REPLACE [+] "x" safe_component ${component} )
   find_library( Lib_${safe_component}
      NAMES rtt_${component}
      PATHS ${DRACO_LIB_DIR} 
            ${DRACO_DIR}/lib
            ${CMAKE_INSTALL_PREFIX}/lib )
   mark_as_advanced(Lib_${safe_component})
endforeach()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(DRACO DEFAULT_MSG DRACO_INCLUDE_DIR Lib_dsxx)

if( DRACO_FOUND )
   get_filename_component( DRACO_LIBRARY_DIR ${Lib_dsxx} PATH )
   set( DRACO_LIBRARY_DIR ${DRACO_LIBRARY_DIR} CACHE PATH 
      "Draco library path." )
   set( DRACO_LIBRARIES "" )
   foreach( component ${components} )
      if( EXISTS ${Lib_${component}} )
         list( APPEND DRACO_LIBRARIES ${Lib_${component}} )
      endif()
   endforeach()
   get_filename_component( DRACO_DIR ${DRACO_LIBRARY_DIR} PATH )
   set( DRACO_DIR ${DRACO_DIR} CACHE PATH "Draco install prefix.")
   message("DRACO_DIR = ${DRACO_DIR}
           " )
# else()

#    message("Draco not found. Try to use cvs checkout.")

#    if( IS_DIRECTORY "/ccs/codes" )
#       set( cvsroot "/ccs/codes/radtran/vendors/cvsroot" )
#    elseif( EXISTS "/usr/projects/draco/yellow" )
#       set( cvsroot "ccscs8:/ccs/codes/radtran/vendors/cvsroot" )
#    else()
#       set( cvsroot "/usr/projects/draco/cvsroot" )
#    endif()

#    message("   CVSROOT = ${cvsroot}")

#    include(ExternalProject)
#    ExternalProject_Add( dracoEP
#       PREFIX         ${PROJECT_BINARY_DIR}/dracoEP
#       CVS_REPOSITORY ${cvsroot}
#       CVS_MODULE     draco
#       UPDATE_COMMAND "cvs -q update -AdP draco"
#       BINARY_DIR     ${PROJECT_BINARY_DIR}/dracoEP )

#    find_path( DRACO_INCLUDE_DIR ds++/config.h
#       ${CMAKE_BINARY_DIR}/dracoEP/include )
#    foreach( component ${components} )
#       message("Looking for Lib_${component} ")
#       find_library( Lib_${component}
#          NAMES rtt_${component}
#          PATHS ${CMAKE_BINARY_DIR}/dracoEP/lib )
#       mark_as_advanced(Lib_${component})
#    endforeach()

#    message("")

endif()


#----------------------------------------------------------------------#
# End 
#----------------------------------------------------------------------#
