#-----------------------------*-cmake-*----------------------------------------#
# file   config/unix-pathf90.cmake
# author Kelly Thompson
# date   2010 Sep 27
# brief  Establish flags for Windows - Intel Visual Fortran
# note   Copyright © 2010 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

# Let anyone who is interested in which FORTRAN compiler we're using
# switch on this macro.
set( CMAKE_Fortran_COMPILER_FLAVOR "PATHSCALE" )

#
# Sanity Checks
# 

if( BUILD_SHARED_LIBS )
  message( FATAL_ERROR "Feature not available - yell at KT." )
endif( BUILD_SHARED_LIBS )

#
# config.h settings
#

# Version information
execute_process( 
  COMMAND ${CMAKE_Fortran_COMPILER} -dumpversion
  OUTPUT_VARIABLE CMAKE_Fortran_COMPILER_VERSION 
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )

#----------------------------------------
# Compiler Flags
#----------------------------------------

# -LNO:simd=0 is a workaround to pathscale 3.1 bug 13958. Les
# Faby, Pathscale 07 November 2007. This sets Loop Nest
# Optimization: vectorization of inner loops to OFF. Without this,
# -O3 produces "invalid read of size 8" warnings under valgrind.  
SET( CMAKE_Fortran_FLAGS 
  "-mcmodel=medium -freeform -OPT:Olimit=0 -fno-second-underscore -LNO:simd=0" )

# -g Provide debug symbols.
# -Wall Issues all compile warnings.
# -O optimized code
# -mtune=opteron optimize binaries for opteron cpus.
SET( CMAKE_Fortran_FLAGS_DEBUG          "-g -Wall -ffortran-bounds-check -DDEBUG" )
SET( CMAKE_Fortran_FLAGS_RELEASE        "-O -mtune=opteron -DNDEBUG" ) # -OPT:Ofast only for 32-bit builds?
SET( CMAKE_Fortran_FLAGS_MINSIZEREL     "${CMAKE_Fortran_FLAGS_RELEASE}" )
SET( CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O -mtune=opteron -g -DDEBUG" ) 

# -mcmodel=medium is not compatible with -fPIC for pathf90.
# Remove -fPIC:
set( CMAKE_SHARED_LIBRARY_Fortran_FLAGS "" )

# ONLY non-debug versions are redistributable.
set( f90_system_dll
  libpathfortran.so
  libopenmp.so  # --> libpthread.so
   )
set( f90_system_lib
  libpathfortran.a
  libopenmp.a
  libpscrt.a
  )

# Static libraries from the /lib directory (useful for
# target_link_library command. 
set( CMAKE_Fortran_compiler_libs "" CACHE INTERNAL 
"Fortran system libraries that are needed by the applications built with Intel Visual Fortran (only optimized versions are redistributable.)" 
)

execute_process(
  COMMAND ${CMAKE_Fortran_COMPILER} --print-search-dirs
  OUTPUT_VARIABLE CMAKE_Fortran_LIB_DIR )
string( REGEX REPLACE ".*libraries[: ]+([/A-z0-9_.]+[^:]).*" "\\1"
  CMAKE_Fortran_LIB_DIR ${CMAKE_Fortran_LIB_DIR} )

foreach( lib ${f90_system_lib} )

  get_filename_component( libwe ${lib} NAME_WE )
  # optimized library
  find_file( CMAKE_Fortran_${libwe}_lib_RELEASE
    NAMES ${lib}
    PATHS "${CMAKE_Fortran_LIB_DIR}"
    )
  mark_as_advanced( CMAKE_Fortran_${libwe}_lib_RELEASE )
  # debug library
  set( CMAKE_Fortran_${libwe}_lib_DEBUG ${CMAKE_Fortran_${libwe}_lib_RELEASE} )
  mark_as_advanced( CMAKE_Fortran_${libwe}_lib_DEBUG )
  set( CMAKE_Fortran_${libwe}_lib_LIBRARY
    optimized
    "${CMAKE_Fortran_${libwe}_lib_RELEASE}"
    debug
    "${CMAKE_Fortran_${libwe}_lib_DEBUG}"
    CACHE INTERNAL "Fortran static system libraries that are needed by the applications built with Intel Visual Fortran (only optimized versions are redistributable.)" FORCE )
  list( APPEND CMAKE_Fortran_compiler_libs ${CMAKE_Fortran_${libwe}_lib_LIBRARY} )

endforeach( lib ${f90_system_lib} )# may need libpthread

# May need to add libgcov here.


if( EXISTS ${CMAKE_Fortran_libopenmp_lib_RELEASE} )

  find_package( Threads )
  list( APPEND CMAKE_Fortran_compiler_libs
    ${CMAKE_THREAD_LIBS_INIT} )

endif( EXISTS ${CMAKE_Fortran_libopenmp_lib_RELEASE} )

#------------------------------------------------------------------------------#
# End config/unix-pathf90.cmake
#------------------------------------------------------------------------------#
