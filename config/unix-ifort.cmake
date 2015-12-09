#-----------------------------*-cmake-*----------------------------------------#
# file   config/unix-ifort.cmake
# author Kelly Thompson
# date   2008 May 30
# brief  Establish flags for Unix - Intel Fortran
# note   Copyright (C) 2010-2015 Los Alamos National Security, LLC.
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
  set( CMAKE_Fortran_COMPILER_VERSION ${CMAKE_Fortran_COMPILER_VERSION} CACHE
    STRING "Fortran compiler version string" FORCE )
  mark_as_advanced( CMAKE_Fortran_COMPILER_VERSION )

  # [KT 2015-07-10] -diag-disable 11060 -- disable warning that is
  #    issued when '-ip' is turned on and a library has no symbols (this
  #    occurs when capsaicin links some trilinos libraries.)
  set( CMAKE_Fortran_FLAGS                "-warn  -fpp -implicitnone -diag-disable 11060" )
  # The cmake configuration for FortanCInterface will only work on
  # Cray CLE if the '-dynamic' option is passed to the linker:
  if( "${CMAKE_Fortran_COMPILER}" MATCHES "ftn" )
    set( CMAKE_Fortran_FLAGS "-dynamic ${CMAKE_Fortran_FLAGS}" )
  endif()
  set( CMAKE_Fortran_FLAGS_DEBUG          "-g -O0 -traceback -ftrapuv -check -DDEBUG" )
  set( CMAKE_Fortran_FLAGS_RELEASE        "-O2 -inline-level=2 -fp-speculation fast -fp-model fast -align array32byte -openmp-report0 -funroll-loops -DNDEBUG" )
  set( CMAKE_Fortran_FLAGS_MINSIZEREL     "${CMAKE_Fortran_FLAGS_RELEASE}" )
  set( CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g -O2 -inline-level=2 -openmp-report0 -inline-level=2 -funroll-loops -DDEBUG" )

  # For older versions of ifort, suppress remarks about vectorization.
  # Note: CMAKE_Fortran_COMPILER_VERSION is not available as of CMake/3.0.0.  We will assume that
  # we are using a Fortran compile1r that is from the same flavor and version as the selected C++
  # compiler.
  if( CMAKE_CXX_COMPILER_VERSION VERSION_LESS 15 )
    set( CMAKE_Fortran_FLAGS_RELEASE        "${CMAKE_Fortran_FLAGS_RELEASE} -vec-report0 -par-report0")
    set( CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO} -vec-report0 -par-report0")
  endif()

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
toggle_compiler_flag( HAVE_MIC    "-mmic"           "Fortran" "")
if( NOT ${SITENAME} STREQUAL "Trinitite" )
  toggle_compiler_flag( ENABLE_SSE  "-mia32 -axSSSE3" "Fortran" "") #sse3, ssse3
endif()
# Use OpenMP_C_FLAGS here because cmake/3.1.1 fails to set
# CMAKE_Fortran_COMPILER_VERSION for FC=mpiifort and FindOpenMP
# chooses the deprecated '-openmp' instead of '-qopenmp'
# Bug reported to Kitware by KT on 2015-01-26
# http://public.kitware.com/Bug/view.php?id=15372
toggle_compiler_flag( OPENMP_FOUND ${OpenMP_C_FLAGS} "Fortran" "" )

# When cross-compiling with '-mmic', rpaths for libraries built from
# Fortran code don't appear to be reported to the icpc linker.  As a
# hack, save these libraries here for manual linking.
if( HAVE_MIC )
  find_library( libifport_loc   NAMES ifport   HINTS ENV LD_LIBRARY_PATH )
  find_library( libifcoremt_loc NAMES ifcoremt HINTS ENV LD_LIBRARY_PATH )
endif()

#------------------------------------------------------------------------------#
# End config/unix-ifort.cmake
#------------------------------------------------------------------------------#
