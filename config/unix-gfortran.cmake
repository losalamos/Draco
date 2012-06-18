#-----------------------------*-cmake-*----------------------------------------#
# file   config/unix-gfortran.cmake
# author Kelly Thompson 
# date   2010 Sep 27
# brief  Establish flags for Windows - Intel Visual Fortran
# note   Copyright (C) 2010-2012 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

# Let anyone who is interested in which FORTRAN compiler we're using
# switch on this macro.
set( CMAKE_Fortran_COMPILER_FLAVOR "GFORTRAN" )

# I know that gfortran 4.1 won't compile our code (maybe 4.2 or 4.3
# will).
message("CMAKE_CXX_COMPILER_VERSION = ${CMAKE_CXX_COMPILER_VERSION}")
if( "${CMAKE_CXX_COMPILER_VERSION}" STRLESS "4.3" )
  message( FATAL_ERROR """
*** Compiler incompatibility:
gfortran < 4.3 will not compile this code.  New versions of gfortran might work but they haven't been tested.  You are trying to use gfortran ${CMAKE_Fortran_COMPILER_VERSION}.
"""
  )
# If we absolutely must compile with an older version of gfortran, the
# following flags may need to be added: "-x f95-cpp-input"

endif()

# General flags:
#
# -ffree-form             Don't assume fixed-form FORTRAN.
# -ffree-line-length-none Allow Fortran lines longer than 132 chars.
# -fimplicit-none         Do not allow implicit typing
# -fPIC                   Produce position independent code
# -cpp                    Enable preprocessing
# -static-libgfortran     On systems that provide libgfortran as a
#                         shared and a static library, this option
#                         forces the use of the static version. If no
#                         shared version of libgfortran was built when
#                         the compiler was configured, this option has
#                         no effect. 

set( CMAKE_Fortran_FLAGS "-ffree-line-length-none -cpp" )

# Debug flags:
#
# -g            Always include debug symbols (strip for distribution).
# -Wall          Enable most warning messages
# -fbounds-check Turn on array bounds checking
# -frange-check  Disallow 1/0=+infinity or range wrap-around.
# -ffpe-trap=invalid  invalid floating point operation, such as "SQRT(-1.0)"
#           =zero     division by zero
#           =overflow overflow in a floating point operation
# -pedantic      Issue warnings for uses of extensions to Fortran 95.
# Consider adding:
# -Wconversion   Warn about implicit conversions between different
#                types.
# -Wimplicit-interface Warn if a procedure is called without an
#                explicit interface.  
#   
SET( CMAKE_Fortran_FLAGS_DEBUG 
  "-g -fbounds-check -frange-check -ffpe-trap=invalid,zero,overflow -fbacktrace -finit-integer=2147483647 -finit-real=NAN -finit-character=127 -DDEBUG" )
# [2011-10-25 KT]: turn off '-Wall -W' for now
# removed pedantic - needed by capsaicin

# Release flags
#
# -03                 Highest supported optimization level.
SET( CMAKE_Fortran_FLAGS_RELEASE "-O3 -ftree-vectorize -funroll-loops -march=k8 -DNDEBUG" )
SET( CMAKE_Fortran_FLAGS_MINSIZEREL "${CMAKE_Fortran_FLAGS_RELEASE}" )
SET( CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g -O3 -ftree-vectorize -funroll-loops -march=k8 -DDEBUG" ) 

# ------------------------------------------------------------
# Find and save compiler libraries.  These may need to be used when
# the main code is C++ that links to Fortran libraries.
# ------------------------------------------------------------

set( f90_system_lib libgfortran.a libgomp.a )
if( ENABLE_OPENMP )
  set( f90_system_lib ${f90_system_lib};libgomp.a )
endif( ENABLE_OPENMP )

# Static libraries from the /lib directory (useful for target_link_library command).
set( CMAKE_Fortran_compiler_libs "" CACHE INTERNAL "Fortran system libraries that are needed by the applications built with Intel Visual Fortran (only optimized versions are redistributable.)" )

execute_process(
  COMMAND ${CMAKE_Fortran_COMPILER} --print-file-name=libgfortran.a
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

# Code Coverage options:

# option( ENABLE_Fortran_CODECOVERAGE
#   "Instrument for Fortran (gfortran) code coverage analysis?" OFF )
# if( ENABLE_Fortran_CODECOVERAGE )
#   find_program( COVERAGE_COMMAND gcov )
#   set( CMAKE_Fortran_FLAGS_DEBUG
#     "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -fprofile-arcs -ftest-coverage" )
#   set( CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS_DEBUG}")
#   set( CMAKE_LDFLAGS               "${CMAKE_LDFLAGS} -fprofile-arcs -ftest-coverage" )
# endif( ENABLE_Fortran_CODECOVERAGE )

# Profiling

# option( ENABLE_Fortran_PROFILING 
#   "Instrument for Fortran profiling with gprof"
#   OFF )
# if( ENABLE_Fortran_PROFILING )
#   set( CMAKE_Fortran_FLAGS   "${CMAKE_Fortran_FLAGS} -pg" )
# endif( ENABLE_Fortran_PROFILING )

#------------------------------------------------------------------------------#
# End config/unix-gfortran.cmake
#------------------------------------------------------------------------------#
