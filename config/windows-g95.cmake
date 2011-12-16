#-----------------------------*-cmake-*----------------------------------------#
# file   config/windows-g95.cmake
# author Kelly Thompson 
# date   2010 May 30
# brief  Establish flags for Windows - Intel Visual Fortran
# note   © Copyright 2010 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

# Let anyone who is interested in which FORTRAN compiler we're using
# switch on this macro.
set( CMAKE_Fortran_COMPILER_FLAVOR "G95" )

# General flags:
#
# -cpp          Preprocess FORTRAN sources.
# -ffree-form   Don't assume fixed-form FORTRAN.
# -ffree-line-length-huge Allow Fortran lines longer than 132 chars.
# -fno-second-underscore Don't allow symbol names in libraries that
#               have 2 trailing underscores.
# -g            Always include debug symbols (strip for distribution).
# -Wno=112 suppress warning #140: "Implicit conversion may cause precision loss."
# -Wno=140 suppress warning #112: "Variable is set but never used."
SET(CMAKE_Fortran_FLAGS 
  "-cpp -ffree-form -ffree-line-length-huge -fno-second-underscore" )

# Debug flags:
#
# -Wall          Enable most warning messages
# -Wextra        Enable warnings not enabled by -Wall
# -ftrace=full   Allow stack tracebacks on abnormal end of program.
# -fbounds-check Turn on array bounds checking
SET(CMAKE_Fortran_FLAGS_DEBUG "-g -Wall -Wextra -fbounds-check -DDEBUG" )

# Release flags
#
# -03                 Highest supported optimization level.
# -funroll-loops      Unroll loops whose number of iterations can be
#                     determined at compile time or upon entry to the
#                     loop. This option makes code larger, and may or
#                     may not make it run faster. 
#-fomit-frame-pointer Don't keep the frame pointer in a register for
#                     functions that don't need one.  
#-msse2
#-mfpmath=sse         Generate floating point arithmetics for SSE
#                     instruction set.
#-mtune=k8            AMD K8 core based CPU with x86_64 instruction
#                     set support.
#      =nocona        Improved version of Intel Pentium4 CPU with
#                     64-bit extensions, MMX, SSE, SSE2 and SSE3
#                     instruction support.
#-march=k8            Implies -mtune=k8
SET(CMAKE_Fortran_FLAGS_RELEASE 
  "-O3 -funroll-loops -fomit-frame-pointer -msse2 -mfpmath=sse -march=k8 -DNDEBUG" )
SET(CMAKE_Fortran_FLAGS_MINSIZEREL "${CMAKE_Fortran_FLAGS_RELEASE}" )
SET(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g -DDEBUG ${CMAKE_Fortran_FLAGS_RELEASE}" ) 

# Save g95 version value
#
execute_process( 
  COMMAND ${CMAKE_Fortran_COMPILER} --version
  OUTPUT_VARIABLE tmp )
string( REGEX REPLACE "^(G95 [(].*)([0-9]+[)]).*" "\\1\\2"
  CMAKE_Fortran_COMPILER_VERSION "${tmp}" )

#
# During discovery of F95 compiler, also discover and make available:
#
# ${CMAKE_Fortran_redist_dll}    
#   - List of Fortran compiler libraries to be installed for release
# ${CMAKE_Fortran_debug_dll}
#   - List of Fortran compiler libraries to be installed for portable developer-only debug version
# ${CMAKE_Fortran_compiler_libs} 
#   - List of Fortran compiler libraries to be used with the target_link_libraries command (C main code that links with Fortran built library.)

# ONLY non-debug versions are redistributable.
set( f90_system_dll
#  none?
   )
set( f90_system_lib
  libf95.a
  # libgcc.a
  )

# Static libraries from the /lib directory (useful for target_link_library command.
set( CMAKE_Fortran_compiler_libs "" CACHE INTERNAL "Fortran system libraries that are needed by the applications built with Intel Visual Fortran (only optimized versions are redistributable.)" )

execute_process(
  COMMAND ${CMAKE_Fortran_COMPILER} --print-file-name=libf95.a
  OUTPUT_VARIABLE CMAKE_Fortran_LIB_DIR )
  get_filename_component( CMAKE_Fortran_LIB_DIR
    ${CMAKE_Fortran_LIB_DIR} PATH )

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

endforeach( lib ${f90_system_lib} )
      
#------------------------------------------------------------------------------#
# End config/windows-ifort.cmake
#------------------------------------------------------------------------------#
