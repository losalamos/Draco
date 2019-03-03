#-----------------------------*-cmake-*----------------------------------------#
# file   FortranChecks/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for FortranChecks.
# note   Copyright (C) 2016-2019, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.9.0)
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
# Build package libraries
# ---------------------------------------------------------------------------- #
add_component_library(
  TARGET       Lib_FC_Derived_Type
  TARGET_DEPS  "Lib_dsxx"
  LIBRARY_NAME FC_Derived_Type
  SOURCES      "${Csources}"
  NOEXPORT )
set_target_properties( Lib_FC_Derived_Type PROPERTIES FOLDER FortranChecks )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #
# Do not install these libraries. They are only used for testing the build
# system.
# install( TARGETS Lib_FC_Derived_Type EXPORT draco-targets DESTINATION lib )

#------------------------------------------------------------------------------#
# Build a library from Fortran sources
#------------------------------------------------------------------------------#

if( BUILD_TESTING )

  if( HAVE_Fortran )
    # For Unix, this will simply run add_subdirectory(f90sub). For Visual Studio
    # or XCode, an ExternalProject will be created that uses an alternate
    # gfortran via Makefiles to generate the library external to the main
    # project.
    include(CMakeAddFortranSubdirectory)

    # CMake does not support storing a list of lists when sending data to a
    # macro.  Because Draco_TPL_INCLUDE_DIRS is a list and we want to stuff it
    # into the list build_system_state, recode Draco_TPL_INCLUDE_DIRS by
    # replacing semicolons with triple underscores.  The list will be
    # reconstructed in the subdirectory's CMakeLists.txt.
    string( REGEX REPLACE ";" "___" tmp
      "${Draco_TPL_INCLUDE_DIRS};${MPI_Fortran_INCLUDE_PATH}")

    # The alternate build system (Makefiles if we are Apple/OSX or Linux/Ninja)
    # will need some of the current build system parameters:
    set( draco_DIR ${Draco_SOURCE_DIR}/config )
    set( build_system_state
#      "-DCMAKE_VERBOSE_MAKEFILE=TRUE"
      "-DDRACO_C4=${DRACO_C4}"
      "-DDRACO_LIBRARY_TYPE=${DRACO_LIBRARY_TYPE}"
      "-DDraco_TPL_INCLUDE_DIRS=${tmp}"
      "-DPYTHON_EXECUTABLE=${PYTHON_EXECUTABLE}"
      "-DRANDOM123_INCLUDE_DIR=${RANDOM123_INCLUDE_DIR}"
      "-Ddraco_DIR=${draco_DIR}"
      )
    if( WIN32 )
      # For Win32 builds, DLL and applications are built in the directory
      # specified by CMAKE_RUNTIME_OUTPUT_DIRECTORY.
      set( build_system_state "${build_system_state};-DCMAKE_RUNTIME_OUTPUT_DIRECTORY=${Draco_BINARY_DIR}/\${CMAKE_BUILD_TYPE}" )
    else()
      set( build_system_state "${build_system_state};-DWITH_CUDA=${WITH_CUDA}" )
    endif()

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
      CMAKE_COMMAND_LINE "${build_system_state}"
      # VERBOSE
      )
  endif()

  add_subdirectory( test )

endif()

# ---------------------------------------------------------------------------- #
# Autodoc
# ---------------------------------------------------------------------------- #
process_autodoc_pages()

# ---------------------------------------------------------------------------- #
# end FortranChecks/CMakeLists.txt
# ---------------------------------------------------------------------------- #
