#-----------------------------*-cmake-*----------------------------------------#
# file   config/unix-pgf90.cmake
# author Kelly Thompson
# date   2011 June 7
# brief  Establish flags for Unix - PGI Fortran
# note   Copyright (C) 2016 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

# Let anyone who is interested in which FORTRAN compiler we're using
# switch on this macro.
set( CMAKE_Fortran_COMPILER_FLAVOR "PGI" )

#
# Compiler Flags
#
if( NOT Fortran_FLAGS_INITIALIZED )
  set( CMAKE_Fortran_COMPILER_VERSION ${CMAKE_Fortran_COMPILER_VERSION} CACHE
       STRING "Fortran compiler version string" FORCE )
  mark_as_advanced( CMAKE_Fortran_COMPILER_VERSION )
  set( Fortran_FLAGS_INITIALIZED "yes" CACHE INTERNAL "using draco settings." )
  set( CMAKE_Fortran_FLAGS                "-Mpreprocess" )
  set( CMAKE_Fortran_FLAGS_DEBUG          "-g -Mbounds -Mchkptr")
  set( CMAKE_Fortran_FLAGS_RELEASE        "-O3")
  set( CMAKE_Fortran_FLAGS_MINSIZEREL     "${CMAKE_Fortran_FLAGS_RELEASE}" )
  set( CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_DEBUG} -O3")
endif()

##---------------------------------------------------------------------------##
# Ensure cache values always match current selection
##---------------------------------------------------------------------------##
set( CMAKE_Fortran_FLAGS                "${CMAKE_Fortran_FLAGS}"                CACHE STRIG "compiler flags" FORCE )
set( CMAKE_Fortran_FLAGS_DEBUG          "${CMAKE_Fortran_FLAGS_DEBUG}"          CACHE STRIG "compiler flags" FORCE )
set( CMAKE_Fortran_FLAGS_RELEASE        "${CMAKE_Fortran_FLAGS_RELEASE}"        CACHE STRIG "compiler flags" FORCE )
set( CMAKE_Fortran_FLAGS_MINSIZEREL     "${CMAKE_Fortran_FLAGS_MINSIZEREL}"     CACHE STRIG "compiler flags" FORCE )
set( CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}" CACHE STRIG "compiler flags" FORCE )

#
# Toggle compiler flags for optional features
#
toggle_compiler_flag( OPENMP_FOUND ${OpenMP_Fortran_FLAGS} "Fortran" "" )

# [2015-01-21 KT] Not sure if we still need the code below...
# ------------------------------------------------------------
# Find and save compiler libraries.  These may need to be used when
# the main code is C++ that links to Fortran libraries.
# ------------------------------------------------------------

# # Order of libraries is important
# set( f90_system_lib libzceh.a libstd.a libC.a )

# # Static libraries from the /lib directory (useful for target_link_library command).
# set( CMAKE_Fortran_compiler_libs "" CACHE INTERNAL
#    "Fortran system libraries that are needed by the applications built with Intel Fortran (only optimized versions are redistributable.)" )

# # Intel Fortran lib directory
# get_filename_component( CMAKE_Fortran_BIN_DIR ${CMAKE_Fortran_COMPILER} PATH )
# string( REPLACE "bin" "lib" CMAKE_Fortran_LIB_DIR ${CMAKE_Fortran_BIN_DIR} )

# # Generate a list of run time libraries.
# foreach( lib ${f90_system_lib} )

#    get_filename_component( libwe ${lib} NAME_WE )
#    # optimized library
#    find_file( CMAKE_Fortran_${libwe}_lib
#       NAMES ${lib}
#       PATHS "${CMAKE_Fortran_LIB_DIR}"
#       HINTS ENV LD_LIBRARY_PATH
#       )
#    mark_as_advanced( CMAKE_Fortran_${libwe}_lib )
#    list( APPEND CMAKE_Fortran_compiler_libs ${CMAKE_Fortran_${libwe}_lib} )

# endforeach()

#------------------------------------------------------------------------------#
# End config/unix-pgf90.cmake
#------------------------------------------------------------------------------#
