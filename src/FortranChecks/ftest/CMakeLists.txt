# ---------------------------------------------------------------------------- #
# Draco - Support libraries for transport solvers.
# Copyright 2012 LANS, LLC.
# ---------------------------------------------------------------------------- #
# $Id: CMakeLists.txt 6136 2011-09-14 00:16:01Z kellyt $
# ---------------------------------------------------------------------------- #

project( FortranChecks_ftest Fortran )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB test_sources *.f90 )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

# BUG: We shouldn't need the special {CMAKE_Fortran_compiler_libs} on
#      the link line, but because of PGI 11-12 broken zc_eh feature,
#      we must for the Fortran compiler to link the nozc_eh C++
#      libraries.  This variable is set in config/unix-pgf90.cmake.

add_scalar_tests( 
   SOURCES "${test_sources}" 
   DEPS    "Lib_FortranChecks;Lib_dsxx;${CMAKE_Fortran_compiler_libs}" )

if( ${CMAKE_Fortran_COMPILER_FLAVOR} STREQUAL "PGI" )
   # CMake has trouble figuring out how PGI works.
   set_target_properties( Ut_FortranChecks_ftest_f90main_exe
      PROPERTIES 
        LINKER_LANGUAGE Fortran )
 endif()

