#----------------------------------------------------------------------#
# file   : platform_checks.cmake
# brief  : Platform Checks for Draco Build System
# note   : Copyright (C) 2011 Los Alamos National Security, LLC.  
#          All rights reserved
# version: $Id$
#----------------------------------------------------------------------#

##---------------------------------------------------------------------------##
## Determine System Type and System Names
##
## Used by ds++ and c4.
##---------------------------------------------------------------------------##
macro( set_draco_uname )
    # Store platform information in config.h
    if( UNIX )
      set( draco_isLinux 1 )
      set( DRACO_UNAME "Linux" )
    elseif( WIN32 )
      set( draco_isWin 1 )
      set( DRACO_UNAME "Windows" )
    elseif( OSF1 )
      set( draco_isOSF1 1 )
      set( DRACO_UNAME "OSF1" )
    elseif( APPLE )
      set( draco_isDarwin 1)
      set( DRACO_UNAME "Darwin" )
    else()
      set( draco_isAIX 1 )
      set( DRACO_UNAME "AIX" )
    endif()
    # Special setup for catamount
    if( ${CMAKE_SYSTEM_NAME} MATCHES "Catamount" )
       set( draco_isLinux_with_aprun 1 )
       set( draco_isCatamount 1 )
       set( DRACO_UNAME "Catamount" )
    endif()
endmacro()

##---------------------------------------------------------------------------##
## Determine if gethostname() is available.
## Determine the value of HOST_NAME_MAX.
## 
## Used by ds++/SystemCall.cc and ds++/path.cc
##---------------------------------------------------------------------------##
macro( query_have_gethostname )
    # Platform checks for gethostname()
    include( CheckIncludeFiles )
    check_include_files( unistd.h HAVE_UNISTD_H )
    check_include_files( limits.h HAVE_LIMITS_H )
    check_include_files( winsock2.h HAVE_WINSOCK2_H )
    check_include_files( direct.h HAVE_DIRECT_H )
    check_include_files( sys/param.h HAVE_SYS_PARAM_H )

    # -------------- Checks for hostname and len(hostname) ---------------- #
    # gethostname()
    include( CheckFunctionExists )
    check_function_exists( gethostname HAVE_GETHOSTNAME )

    # HOST_NAME_MAX
    include( CheckSymbolExists )
    unset( hlist )
    if( HAVE_UNISTD_H )
       list( APPEND hlist unistd.h )
    endif()
    if( HAVE_WINSOCK2_H )
       list( APPEND hlist winsock2.h )
    endif()
    if( HAVE_LIMITS_H )
       list( APPEND hlist limits.h )
    endif()
    check_symbol_exists( HOST_NAME_MAX "${hlist}" HAVE_HOST_NAME_MAX )
    if( NOT HAVE_HOST_NAME_MAX )
       unset( HAVE_GETHOSTNAME )
    endif()

    check_symbol_exists( _POSIX_HOST_NAME_MAX "posix1_lim.h" HAVE_POSIX_HOST_NAME_MAX )

    # HOST_NAME_MAX
    check_symbol_exists( MAXHOSTNAMELEN "sys/param.h" HAVE_MAXHOSTNAMELEN )
    if( NOT HAVE_MAXHOSTNAMELEN )
       unset( HAVE_MAXHOSTNAMELEN )
    endif()
    
    # MAX_COMPUTERNAME_LENGTH
    # check_symbol_exists( MAX_COMPUTERNAME_LENGTH "windows.h" HAVE_MAX_COMPUTERNAME_LENGTH )
    # if( NOT HAVE_MAX_COMPUTERNAME_LENGTH )
    #    unset( HAVE_MAX_COMPUTERNAME_LENGTH )
    # endif()
endmacro()

##---------------------------------------------------------------------------##
## Determine if gethostname() is available.
## Determine the value of HOST_NAME_MAX.
## 
## Used by ds++/SystemCall.cc and ds++/path.cc
##---------------------------------------------------------------------------##
macro( query_have_maxpathlen )
    # MAXPATHLEN
    unset( hlist )
    if( HAVE_UNISTD_H )
       list( APPEND hlist unistd.h )
    endif()
    if( HAVE_LIMITS_H )
       list( APPEND hlist limits.h )
    endif()
    if( HAVE_SYS_PARAM_H )
       list( APPEND hlist sys/param.h )
    endif()
    check_symbol_exists( MAXPATHLEN "${hlist}" HAVE_MAXPATHLEN )
    if( NOT HAVE_MAXPATHLEN )
        unset( HAVE_MAXPATHLEN )
    endif()
endmacro()

##---------------------------------------------------------------------------##
## Determine Word Types
##
## Query sizes of PODTs.
##---------------------------------------------------------------------------##

macro( determine_word_types )
   include(CheckTypeSize)
   check_type_size( "int"       SIZEOF_INT )
   check_type_size( "long"      SIZEOF_LONG )
   check_type_size( "long long" SIZEOF_LONG_LONG )
   
   check_type_size( "float"       SIZEOF_FLOAT )
   check_type_size( "double"      SIZEOF_DOUBLE )
   check_type_size( "long double" SIZEOF_LONG_DOUBLE )
endmacro()

##---------------------------------------------------------------------------##
## Check 8-byte int type
##
## For some systems, provide special compile flags to support 8-byte integers
##---------------------------------------------------------------------------##

macro(check_eight_byte_int_type)
   if( "${SIZEOF_INT}notset" STREQUAL "notset" )
      determine_word_types()
   endif()

   if( "${SIZEOF_INT}" STREQUAL "8" )
      message( "Checking for 8-byte integer type... int - no mods needed." )
   elseif( "${SIZEOF_LONG}" STREQUAL "8" )
      message( "Checking for 8-byte integer type... long - no mods needed." )
   else()
      message( FATAL_ERROR "need to patch up this part of the build system." )
# See ac_dracotests.m4
#      if( ${CMAKE_SYSTEM} MATCHES AIX )
#         if with_cxx = ibm and vendor_mpi and strict ansi, then add
#      -qlanglvl=extended -qlonglong
#         fi
#      endif()
   endif()
endmacro()

##---------------------------------------------------------------------------##
## Wedgehog types
##---------------------------------------------------------------------------##

macro( wedgehog_types )
   if( "${HOST_INT_SIZE}" STREQUAL "${SIZEOF_INT}" )
      set( HOST_INT "int" )
   elseif( "${HOST_INT_SIZE}" STREQUAL "${SIZEOF_LONG}" )
      set( HOST_INT "long" )
   elseif( "${HOST_INT_SIZE}" STREQUAL "${SIZEOF_LONG_LONG}" )
      set( HOST_INT "long long" )
   endif()

   if( "${EIGHT_BYTE_INT_SIZE}" STREQUAL "${SIZEOF_INT}" )
      set( EIGHT_BYTE_INT "int" )
   elseif( "${EIGHT_BYTE_INT_SIZE}" STREQUAL "${SIZEOF_LONG}" )
      set( EIGHT_BYTE_INT "long" )
   elseif( "${EIGHT_BYTE_INT_SIZE}" STREQUAL "${SIZEOF_LONG_LONG}" )
      set( EIGHT_BYTE_INT "long long" )
   endif()

   if( "${HOST_FLOAT_SIZE}" STREQUAL "${SIZEOF_FLOAT}" )
      set( HOST_FLOAT "float" )
   elseif( "${HOST_FLOAT_SIZE}" STREQUAL "${SIZEOF_DOUBLE}" )
      set( HOST_FLOAT "double" )
   endif()
endmacro()

##---------------------------------------------------------------------------##



