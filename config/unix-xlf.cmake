#--------------------------------------------*-cmake-*---------------------------------------------#
# file   config/unix-xlf.cmake
# author Gabriel Rockefeller, Kelly Thompson <kgt@lanl.gov>
# date   2012 Nov 1
# brief  Establish flags for Unix - IBM XL Fortran
# note   Copyright (C) 2016-2020 Triad National Security, LLC.
#        All rights reserved.
#--------------------------------------------------------------------------------------------------#

# Let anyone who is interested in which FORTRAN compiler we're using switch on
# this macro.
set( CMAKE_Fortran_COMPILER_FLAVOR "XL" )

# -qsuppress=
#   1501-210 (W) command option t contains an incorrect subargument
#            This is related to '-pthread' showing up in the flags (a cmake bug)

if( NOT Fortran_FLAGS_INITIALIZED )
   set( Fortran_FLAGS_INITIALIZED "yes" CACHE INTERNAL "using draco settings." )
   set( CMAKE_Fortran_FLAGS
     "-qlanglvl=2008std -qflag=i:w -qsuppress=1501-210" )
   if( ${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "ppc64le")
     string( APPEND CMAKE_Fortran_FLAGS " -qarch=pwr9 -qtune=pwr9" )
   endif()

   set( CMAKE_Fortran_FLAGS_DEBUG "-O0" )
   set( CMAKE_Fortran_FLAGS_RELWITHDEBINFO
     "-O3 -qstrict=nans:operationprecision" )
   set( CMAKE_Fortran_FLAGS_RELEASE
     "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}" )
   set( CMAKE_Fortran_FLAGS_MINSIZEREL "${CMAKE_Fortran_FLAGS_RELEASE}" )
endif()

##---------------------------------------------------------------------------##
# Ensure cache values always match current selection
##---------------------------------------------------------------------------##
set( CMAKE_Fortran_FLAGS                "${CMAKE_Fortran_FLAGS}"
  CACHE STRING "compiler flags" FORCE )
set( CMAKE_Fortran_FLAGS_DEBUG          "${CMAKE_Fortran_FLAGS_DEBUG}"
  CACHE STRING "compiler flags" FORCE )
set( CMAKE_Fortran_FLAGS_RELEASE        "${CMAKE_Fortran_FLAGS_RELEASE}"
  CACHE STRING "compiler flags" FORCE )
set( CMAKE_Fortran_FLAGS_MINSIZEREL     "${CMAKE_Fortran_FLAGS_MINSIZEREL}"
  CACHE STRING "compiler flags" FORCE )
set( CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}"
  CACHE STRING "compiler flags" FORCE )

#--------------------------------------------------------------------------------------------------#
# End config/unix-xlf.cmake
#--------------------------------------------------------------------------------------------------#
