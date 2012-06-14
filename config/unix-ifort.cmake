#-----------------------------*-cmake-*----------------------------------------#
# file   config/unix-ifort.cmake
# author Kelly Thompson 
# date   2008 May 30
# brief  Establish flags for Unix - Intel Fortran
# note   Copyright (C) 2010-2012 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

# Let anyone who is interested in which FORTRAN compiler we're using
# switch on this macro.
set( CMAKE_Fortran_COMPILER_FLAVOR "IFORT" )

# General flags:
#
# -W1           enable all warnings.
# -free         free format FORTRAN
# -fpp          this file.
# -fp-model strict - Enables value-safe optimizations on floating-point data 
#              and rounds intermediate results to source-defined precision.
# -prec-div    Attempts to use slower but more accurate implementation of 
#              floating-point divide. 
# -static      Link against static Fortran runtime libraries.
# -fPIC        Generate position independent code
# -implicitnone Do not allow automatic variable types.
# -openmp      Enable OpenMP parallelization
# -parallel    Automatic parallelziation.
set( CMAKE_Fortran_FLAGS 
  "-warn  -fpp -fPIC -implicitnone"
  ) #  -fp-model strict -prec-div -static
if( USE_OPENMP )
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -openmp -parallel" )
endif()

# Debug flags:
#
# -g            Always include debug symbols (strip for distribution).
# -O0           Disable optimizations
# -traceback    Provide extra traceback information.
# -ftrapuv      Initialize stack local variables to an unusual value.
# -check bounds Turn on array bounds checking
# -check uninit Check for uninitialized variables
# -check        Turn on all checks
# -check pointers run-time checking for disassociated or uninitialized pointers.

set( CMAKE_Fortran_FLAGS_DEBUG 
  "-g -O0 -traceback -ftrapuv -check -DDEBUG" )

# Release flags
#
# -02             optimization level.
# -inline-level=2 specifies the level of inline function expansion.
# -funroll-loops  unroll user loops based on the default optimziation.
# -mtune=core2
# -vec-report0    silence reporting
# -openmp-report0 silence reporting
# -par-report0    silence reporting
SET( CMAKE_Fortran_FLAGS_RELEASE 
  "-O2 -inline-level=2 -vec-report0 -openmp-report0 -par-report0 -funroll-loops -DNDEBUG" )
SET( CMAKE_Fortran_FLAGS_MINSIZEREL "${CMAKE_Fortran_FLAGS_RELEASE}" )
SET( CMAKE_Fortran_FLAGS_RELWITHDEBINFO 
  "-g -O2 -inline-level=2 -vec-report0 -openmp-report0 -par-report0 -inline-level=2 -funroll-loops -DDEBUG" ) 

if( ENABLE_SSE )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mia32 -axSSSE3" ) # sse3, ssse3
endif( ENABLE_SSE )


# Save ifort version value
# execute_process( 
#   COMMAND ${CMAKE_Fortran_COMPILER} --version
#   OUTPUT_VARIABLE tmp )
# string( REGEX REPLACE ".*([0-9]+[.][0-9][ ][0-9]+).*" 
#   "\\1" CMAKE_Fortran_COMPILER_VERSION "${tmp}" 
#   )   

# ------------------------------------------------------------
# Find and save compiler libraries.  These may need to be used when
# the main code is C++ that links to Fortran libraries.
# ------------------------------------------------------------

# Order of libraries is important
set( f90_system_lib 
   libifport.a libifcore.a libirc.a libsvml.a libimf.a )
if( USE_OPENMP )
  set( f90_system_lib ${f90_system_lib};libiomp5.a )
endif( USE_OPENMP )

# Static libraries from the /lib directory (useful for target_link_library command).
set( CMAKE_Fortran_compiler_libs "" CACHE INTERNAL
   "Fortran system libraries that are needed by the applications built with Intel Fortran (only optimized versions are redistributable.)" )

# Intel Fortran lib directory
get_filename_component( CMAKE_Fortran_BIN_DIR ${CMAKE_Fortran_COMPILER} PATH )
string( REPLACE "bin" "lib" CMAKE_Fortran_LIB_DIR ${CMAKE_Fortran_BIN_DIR} )

# Generate a list of run time libraries.
foreach( lib ${f90_system_lib} )

   get_filename_component( libwe ${lib} NAME_WE )
   # optimized library
   find_file( CMAKE_Fortran_${libwe}_lib_RELEASE
      NAMES ${lib}
      PATHS "${CMAKE_Fortran_LIB_DIR}"
      PATH_SUFFIXES "intel64"
      )
   mark_as_advanced( CMAKE_Fortran_${libwe}_lib_RELEASE )
   # debug library
   set( CMAKE_Fortran_${libwe}_lib_DEBUG ${CMAKE_Fortran_${libwe}_lib_RELEASE} )
   mark_as_advanced( CMAKE_Fortran_${libwe}_lib_DEBUG )
   set( CMAKE_Fortran_${libwe}_lib_LIBRARY
      optimized "${CMAKE_Fortran_${libwe}_lib_RELEASE}"
      debug     "${CMAKE_Fortran_${libwe}_lib_DEBUG}"
      CACHE INTERNAL "Fortran static system libraries that are needed by the applications built with Intel Visual Fortran (only optimized versions are redistributable.)" FORCE )
   list( APPEND CMAKE_Fortran_compiler_libs ${CMAKE_Fortran_${libwe}_lib_LIBRARY} )
   
endforeach()

# # Mixed Language Linking:
# # When linking C++ code that calls F90 libraries, we also need -ldl -lpthread

# if( ${CMAKE_CXX_COMPILER} MATCHES "g[+][+]" )

#    # Threads
#    if( NOT CMAKE_THREAD_LIBS_INIT )
#       find_package( Threads )
#    endif()
#    if( CMAKE_THREAD_LIBS_INIT ) 
#       # set( GCC_LIBRARIES ${GCC_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} )
#       set( CMAKE_Fortran_compiler_libs 
#          ${CMAKE_Fortran_compiler_libs} ${CMAKE_THREAD_LIBS_INIT} )
#    endif()
   
#    # libdl.a
#    execute_process( 
#       COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=libdl.so
#       TIMEOUT 5
#       RESULT_VARIABLE tmp
#       OUTPUT_VARIABLE libdl_so_loc
#       ERROR_VARIABLE err
#       OUTPUT_STRIP_TRAILING_WHITESPACE
#       )
#    get_filename_component( libdl_so_loc ${libdl_so_loc} ABSOLUTE )
#    # set( GCC_LIBRARIES ${GCC_LIBRARIES} ${libdl_a_loc} )
#    set( CMAKE_Fortran_compiler_libs ${CMAKE_Fortran_compiler_libs} ${libdl_so_loc} )
   
# endif()

#------------------------------------------------------------------------------#
# End config/unix-ifort.cmake
#------------------------------------------------------------------------------#
