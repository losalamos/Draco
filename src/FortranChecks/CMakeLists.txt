#-----------------------------*-cmake-*----------------------------------------#
# file   config/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for FortranChecks.
# note   Copyright (C) 2010-2012, Los Alamos National Security
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#
cmake_minimum_required( VERSION 2.8 )
project( FortranChecks CXX Fortran )

# When using CXX/Fortran cross language linking and function calling
# the CMake module FortranCInterface.cmake can be used to deal with
# name mangling issues.  However, Draco has chosen to use
# ISO_C_BINDING instead and this makes manual name demangling
# unnecessary.

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB f90sources *.f90 )
file( GLOB Csources *.cc )

# ---------------------------------------------------------------------------- #
# Build package libraries
# ---------------------------------------------------------------------------- #

add_component_library( 
   TARGET       Lib_FortranChecks
   LIBRARY_NAME FortranCheck
   SOURCES      "${f90sources}" 
   LINK_LANGUAGE "Fortran" )

add_component_library( 
   TARGET       Lib_FC_Derived_Type 
   LIBRARY_NAME FC_Derived_Type
   SOURCES      "${Csources}" )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_FortranChecks Lib_FC_Derived_Type
   EXPORT draco-targets DESTINATION lib )

# ---------------------------------------------------------------------------- #
# Unit tests
# ---------------------------------------------------------------------------- #

if( BUILD_TESTING )
   add_subdirectory( test )
   add_subdirectory( ftest )
endif()   

# ---------------------------------------------------------------------------- #
# Autodoc
# ---------------------------------------------------------------------------- #

process_autodoc_pages()

