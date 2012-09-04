# ---------------------------------------------------------------------------- #
# Draco - Support libraries for transport solvers.
# Copyright 2012 LANS, LLC.
# ---------------------------------------------------------------------------- #
# $Id$
# ---------------------------------------------------------------------------- #

project( FortranChecks_ftest Fortran )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB test_sources *.f90 )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #
set( test_deps "Lib_FortranChecks;Lib_dsxx;Lib_FC_Derived_Type" )

# BUG: We shouldn't need the special {CMAKE_Fortran_compiler_libs} on
#      the link line, but because of PGI 11-12 broken zc_eh feature,
#      we must for the Fortran compiler to link the nozc_eh C++
#      libraries.  This variable is set in config/unix-pgf90.cmake.
if( "${CMAKE_Fortran_COMPILER_FLAVOR}" STREQUAL "PGI" )
   set( test_deps "${test_deps};${CMAKE_Fortran_compiler_libs}" )
endif()

# Gfortran prior to 4.7 does not implement iso_c_binding correctly:
if( "${CMAKE_Fortran_COMPILER_FLAVOR}" STREQUAL "GFORTRAN" AND
      "${CMAKE_Fortran_COMPILER_VERSION}" STRLESS "4.7")
   # do not build these tests.
   message( "WARNING: Using gfortran < 4.7.  FortranChecks/ftest disabled.")
else()
   add_scalar_tests( 
      SOURCES "${test_sources}" 
      DEPS    "${test_deps}" )
endif()

# Static library builds and FC=pgfortran builds need help determining
# the linker language.
foreach( ut f90main derived_types )
   set_target_properties( Ut_FortranChecks_ftest_${ut}_exe
      PROPERTIES LINKER_LANGUAGE Fortran )
endforeach()
