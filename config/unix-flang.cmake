#-----------------------------*-cmake-*----------------------------------------#
# file   config/unix-flang.cmake
# author Kelly Thompson
# date   Sunday, Apr 29, 2018, 19:56 pm
# brief  Establish flags for Unix/Linux - Gnu Fortran
# note   Copyright (C) 2018-2019 Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

if( NOT Fortran_FLAGS_INITIALIZED )
  mark_as_advanced( CMAKE_Fortran_COMPILER_VERSION )
endif()

#
# Compiler Flags
#

if( NOT Fortran_FLAGS_INITIALIZED )
   set( Fortran_FLAGS_INITIALIZED "yes" CACHE INTERNAL "using draco settings." )

   set( CMAKE_Fortran_FLAGS "-cpp" )
   set( CMAKE_Fortran_FLAGS_DEBUG "-g -gdwarf-3 -DDEBUG" )
   set( CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG" )
   set( CMAKE_Fortran_FLAGS_MINSIZEREL "${CMAKE_Fortran_FLAGS_RELEASE}" )
   set( CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g -gdwarf-3 ${CMAKE_Fortran_FLAGS_RELEASE}")

endif()

##---------------------------------------------------------------------------##
# Ensure cache values always match current selection
##---------------------------------------------------------------------------##
set( CMAKE_Fortran_FLAGS                "${CMAKE_Fortran_FLAGS}"                CACHE STRING "compiler flags" FORCE )
set( CMAKE_Fortran_FLAGS_DEBUG          "${CMAKE_Fortran_FLAGS_DEBUG}"          CACHE STRING "compiler flags" FORCE )
set( CMAKE_Fortran_FLAGS_RELEASE        "${CMAKE_Fortran_FLAGS_RELEASE}"        CACHE STRING "compiler flags" FORCE )
set( CMAKE_Fortran_FLAGS_MINSIZEREL     "${CMAKE_Fortran_FLAGS_MINSIZEREL}"     CACHE STRING "compiler flags" FORCE )
set( CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}" CACHE STRING "compiler flags" FORCE )

#
# Toggle compiler flags for optional features
#
if( OpenMP_Fortran_FLAGS )
  toggle_compiler_flag( OPENMP_FOUND ${OpenMP_Fortran_FLAGS} "Fortran" "" )
endif()

#------------------------------------------------------------------------------#
# End config/unix-flang.cmake
#------------------------------------------------------------------------------#
