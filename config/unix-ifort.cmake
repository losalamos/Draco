#-----------------------------*-cmake-*----------------------------------------#
# file   config/unix-ifort.cmake
# author Kelly Thompson 
# date   2008 May 30
# brief  Establish flags for Unix - Intel Fortran
# note   Copyright (C) 2010-2013 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

# Let anyone who is interested in which FORTRAN compiler we're using
# switch on this macro.
set( CMAKE_Fortran_COMPILER_FLAVOR "IFORT" )

#
# Compiler flags:
#
if( NOT Fortran_FLAGS_INITIALIZED )
   set( Fortran_FLAGS_INITIALIZED "yes" CACHE INTERNAL "using draco settings." )

   set( CMAKE_Fortran_FLAGS                "-warn  -fpp -implicitnone" ) 
   set( CMAKE_Fortran_FLAGS_DEBUG          "-g -O0 -traceback -ftrapuv -check -DDEBUG" )
   set( CMAKE_Fortran_FLAGS_RELEASE        "-O2 -inline-level=2 -fp-speculation fast -fp-model fast -align array32byte -vec-report0 -openmp-report0 -par-report0 -funroll-loops -DNDEBUG" )
   set( CMAKE_Fortran_FLAGS_MINSIZEREL     "${CMAKE_Fortran_FLAGS_RELEASE}" )
   set( CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g -O2 -inline-level=2 -vec-report0 -openmp-report0 -par-report0 -inline-level=2 -funroll-loops -DDEBUG" ) 

endif()

##---------------------------------------------------------------------------##
# Ensure cache values always match current selection
##---------------------------------------------------------------------------##
set( CMAKE_Fortran_FLAGS                "${CMAKE_Fortran_FLAGS}"                CACHE STRING "compiler flags" FORCE )
set( CMAKE_Fortran_FLAGS_DEBUG          "${CMAKE_Fortran_FLAGS_DEBUG}"          CACHE STRING "compiler flags" FORCE ) 
set( CMAKE_Fortran_FLAGS_RELEASE        "${CMAKE_Fortran_FLAGS_RELEASE}"        CACHE STRING "compiler flags" FORCE )
set( CMAKE_Fortran_FLAGS_MINSIZEREL     "${CMAKE_Fortran_FLAGS_MINSIZEREL}"     CACHE STRING "compiler flags" FORCE )
set( CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}" CACHE STRING "compiler flags" FORCE )

# Optional compiler flags
toggle_compiler_flag( USE_OPENMP  "-openmp" "Fortran" "")
toggle_compiler_flag( ENABLE_SSE  "-mia32 -axSSSE3" "Fortran" "") # sse3, ssse3

#------------------------------------------------------------------------------#
# End config/unix-ifort.cmake
#------------------------------------------------------------------------------#
