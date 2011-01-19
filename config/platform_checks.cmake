#----------------------------------------------------------------------#
# file   : platform_checks.cmake
# author : Kelly Thompson <kgt@lanl.gov>
# date   : 2011 Jan 19  8:58 am
# brief  : Platform Checks for Draco Build System
# note   : Copyright (C) 2011 Los Alamos National Security, LLC.  
#          All rights reserved
# version: $Id$
#----------------------------------------------------------------------#


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
