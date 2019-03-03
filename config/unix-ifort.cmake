#-----------------------------*-cmake-*----------------------------------------#
# file   config/unix-ifort.cmake
# author Kelly Thompson
# date   2008 May 30
# brief  Establish flags for Unix - Intel Fortran
# note   Copyright (C) 2016-2019 Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

#
# Compiler flags:
#
if( NOT Fortran_FLAGS_INITIALIZED )
  set( Fortran_FLAGS_INITIALIZED "yes" CACHE INTERNAL "using draco settings." )
  set( CMAKE_Fortran_COMPILER_VERSION ${CMAKE_Fortran_COMPILER_VERSION} CACHE
    STRING "Fortran compiler version string" FORCE )
  mark_as_advanced( CMAKE_Fortran_COMPILER_VERSION )

  # [KT 2015-07-10] -diag-disable 11060 -- disable warning that is issued when
  #    '-ip' is turned on and a library has no symbols (this occurs when
  #    capsaicin links some trilinos libraries.)
  # [KT 2016-11-16] -diag-disable 11021 -- disable warning that is issued when
  #    '-ip' is turned on and a library has unresolved symbols (this occurs when
  #    capsaicin links to openmpi/1.10.3 on snow/fire/ice).
  #    Ref: https://github.com/open-mpi/ompi/issues/251
  # [KT 2018-03-14] '-assume nostd_mod_proc_name' --  discussion with G.
  #    Rockefeller and S. Nolen aobut ifort's non-standard name mangling for
  #    module procedures. Not sure if we need this yet.

  set( CMAKE_Fortran_FLAGS
    "-warn  -fpp -implicitnone -diag-disable 11060" )
  set( CMAKE_Fortran_FLAGS_DEBUG
    "-g -O0 -traceback -ftrapuv -check -DDEBUG" )
  set( CMAKE_Fortran_FLAGS_RELEASE
    "-O2 -inline-level=2 -fp-speculation fast -fp-model fast" )
  string(APPEND CMAKE_Fortran_FLAGS_RELEASE
    " -align array32byte -funroll-loops -diag-disable 11021 -DNDEBUG" )
  set( CMAKE_Fortran_FLAGS_MINSIZEREL
    "${CMAKE_Fortran_FLAGS_RELEASE}" )
  set( CMAKE_Fortran_FLAGS_RELWITHDEBINFO
    "-g -O2 -inline-level=2 -funroll-loops -DDEBUG" )

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

# Optional compiler flags
if( NOT ${SITENAME} STREQUAL "Trinitite" )
  toggle_compiler_flag( ENABLE_SSE  "-mia32 -axSSSE3" "Fortran" "") #sse3, ssse3
endif()
toggle_compiler_flag( OPENMP_FOUND ${OpenMP_Fortran_FLAGS} "Fortran" "" )

#------------------------------------------------------------------------------#
# End config/unix-ifort.cmake
#------------------------------------------------------------------------------#
