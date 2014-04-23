#-----------------------------*-cmake-*----------------------------------------#
# file   config/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2013 June 24
# brief  Generate build unit tests for fortran/C binding.
# note   Copyright (C) 2013-2014, Los Alamos National Security, LLC.
#        All rights reserved.
# ---------------------------------------------------------------------------- #
# $Id$
# ---------------------------------------------------------------------------- #
cmake_minimum_required(VERSION 2.6)
project( FortranChecks_ftest Fortran )

# ---------------------------------------------------------------------------- #
# Setup - This folder is parsed independent of parent folders, so we need 
# some basic setup.
# ---------------------------------------------------------------------------- #
# Build system configuration files are located here.
if( NOT EXISTS ${PROJECT_SOURCE_DIR}/../../../config )
   message( FATAL_ERROR "can't find ${PROJECT_SOURCE_DIR}/../../../config" )
else()
   set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../../../config )
endif()

# Provide helper functions used by component CMakeLists.txt files.
include( component_macros )
# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #
file( GLOB test_sources *.f90 )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #
set( test_deps "Lib_FC_Derived_Type;Lib_dsxx" )

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
      
# : cannot find -lLib_FortranChecks
# : cannot find -lLib_dsxx
# : cannot find -lLib_FC_Derived_Type
      

   # Static library builds and FC=pgfortran builds need help determining
   # the linker language.
   foreach( ut f90main derived_types )
      set_target_properties( Ut_FortranChecks_ftest_${ut}_exe
         PROPERTIES LINKER_LANGUAGE Fortran )
   endforeach()
endif()
