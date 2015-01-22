#-----------------------------*-cmake-*----------------------------------------#
# file   FortranChecks/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for FortranChecks.
# note   Copyright (C) 2010-2014, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.0.0)
project( FortranChecks CXX )

# When using CXX/Fortran cross language linking and function calling the CMake
# module FortranCInterface.cmake can be used to deal with name mangling issues.
# However, Draco has chosen to use ISO_C_BINDING instead and this makes manual
# name demangling unnecessary.

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #
file( GLOB Csources *.cc )

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories(
   ${dsxx_BINARY_DIR} ) # ds++/config.h

# ---------------------------------------------------------------------------- #
# Build package libraries
# ---------------------------------------------------------------------------- #

add_component_library(
   TARGET       Lib_FC_Derived_Type
   LIBRARY_NAME FC_Derived_Type
   SOURCES      "${Csources}"
   NOEXPORT )

set_target_properties( Lib_FC_Derived_Type
   PROPERTIES FOLDER FortranChecks )

# For Unix, this will simply run add_subdirectory(f90sub). For MVSE or XCode, an
# ExternalProject will be created that uses an alternate gfortran via Makefiles
# to generate the library external to the main project.
include(CMakeAddFortranSubdirectory)

# CMake will put some run scripts in $build_dir/f90sub Another directory for
# the build is generated at $build_dir/FortranChecks_f90sub_build-prefix
cmake_add_fortran_subdirectory(
  f90sub                       # directory name
  PROJECT FortranChecks_f90sub # project name in toplevel CMakeLists.txt
                               # Creates target named 'FortranChecks_f90sub_build'
  ARCHIVE_DIR  f90sub          # .lib location relative to root binary tree
  RUNTIME_DIR  f90sub          # .dll location relative to root binary tree
  LIBRARIES    rtt_FC_f90sub   # Sub project will create a library with this name.
                               # --> librtt_FC_f90.dll.
  TARGET_NAMES Lib_FC_f90sub
  DEPENDS      "Lib_FC_Derived_Type;Lib_dsxx"
  NO_EXTERNAL_INSTALL
# VERBOSE
)

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #
# No need to install this library as it is not used by anyone else.  This
# component is only used to test Fortran/C binding.
# install( TARGETS Lib_FC_Derived_Type EXPORT draco-targets DESTINATION lib )

# ---------------------------------------------------------------------------- #
# Unit tests
# ---------------------------------------------------------------------------- #
if( BUILD_TESTING )
   add_subdirectory( test )
endif()

# ---------------------------------------------------------------------------- #
# Autodoc
# ---------------------------------------------------------------------------- #
process_autodoc_pages()
