#-----------------------------*-cmake-*----------------------------------------#
# file   config/unix-crayftn.cmake
# author Kelly Thompson
# date   2008 May 30
# brief  Establish flags for Unix - Intel Fortran
# note   Copyright (C) 2016-2020 Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

#
# Compiler flags:
#
if( NOT Fortran_FLAGS_INITIALIZED )
  set( Fortran_FLAGS_INITIALIZED "yes" CACHE INTERNAL "using draco settings." )

  # Find and cache the compiler version (2015-09-28 CMake-3.3.1 misses this).
  execute_process( COMMAND ${CMAKE_Fortran_COMPILER} -V
    ERROR_VARIABLE ftn_version_output
    OUTPUT_STRIP_TRAILING_WHITESPACE )
  string( REGEX REPLACE ".*Version ([0-9]+)[.]([0-9]+)[.]([0-9]+).*" "\\1.\\2"
    CMAKE_Fortran_COMPILER_VERSION "${ftn_version_output}" )
  set( CMAKE_Fortran_COMPILER_VERSION ${CMAKE_Fortran_COMPILER_VERSION} CACHE
    STRING "Fortran compiler version string" FORCE )
  mark_as_advanced( CMAKE_Fortran_COMPILER_VERSION )

  set( CMAKE_Fortran_FLAGS                "" )
  set( CMAKE_Fortran_FLAGS_DEBUG          "-g -O0 -DDEBUG" )
  set( CMAKE_Fortran_FLAGS_RELEASE        "-O3 -DNDEBUG" )
  set( CMAKE_Fortran_FLAGS_MINSIZEREL     "${CMAKE_Fortran_FLAGS_RELEASE}" )
  set( CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g -O3 -DDEBUG" )

endif()

##---------------------------------------------------------------------------##
# Ensure cache values always match current selection
##---------------------------------------------------------------------------##
set( CMAKE_Fortran_FLAGS                "${CMAKE_Fortran_FLAGS}"                CACHE STRING "compiler flags" FORCE )
set( CMAKE_Fortran_FLAGS_DEBUG          "${CMAKE_Fortran_FLAGS_DEBUG}"          CACHE STRING "compiler flags" FORCE )
set( CMAKE_Fortran_FLAGS_RELEASE        "${CMAKE_Fortran_FLAGS_RELEASE}"        CACHE STRING "compiler flags" FORCE )
set( CMAKE_Fortran_FLAGS_MINSIZEREL     "${CMAKE_Fortran_FLAGS_MINSIZEREL}"     CACHE STRING "compiler flags" FORCE )
set( CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}" CACHE STRING "compiler flags" FORCE )

toggle_compiler_flag( OPENMP_FOUND ${OpenMP_Fortran_FLAGS} "Fortran" "" )

# -craype-verbose
# -hnegmsgs # show pos/neg messages about optimizations.
# -hlist=m  # creates annotated listing (loopmark).
#


#------------------------------------------------------------------------------#
# End config/unix-crayftn.cmake
#------------------------------------------------------------------------------#
