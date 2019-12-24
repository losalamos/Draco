#-----------------------------*-cmake-*----------------------------------------#
# file   quadrature/CMakeLists.txt
# brief  Instructions for building quadrature Makefiles.
# note   Copyright (C) 2010-2016-2019 Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.9.0)
project( quadrature CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB sources *.cc )
file( GLOB headers *.hh )
file( GLOB f90sources *.f90 )

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library(
   TARGET       Lib_quadrature
   TARGET_DEPS  "Lib_special_functions;Lib_parser;Lib_mesh_element"
   LIBRARY_NAME ${PROJECT_NAME}
   SOURCES      "${sources}"
   HEADERS      "${headers};${f90sources}"
)

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_quadrature EXPORT draco-targets DESTINATION
  ${DBSCFGDIR}lib )
install( FILES ${headers} ${f90sources}
  DESTINATION ${DBSCFGDIR}include/quadrature )
if( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  install(FILES $<TARGET_PDB_FILE:Lib_quadrature> DESTINATION ${DBSCFGDIR}lib
    OPTIONAL)
endif()

# ---------------------------------------------------------------------------- #
# Unit tests
# ---------------------------------------------------------------------------- #

if( BUILD_TESTING )
   add_subdirectory( test )

  if( HAVE_Fortran )
    # For Unix, this will simply run add_subdirectory(ftest). For MVSE or XCode,
    # an ExternalProject will be created that uses an alternate gfortran via
    # Makefiles to generate the library external to the main project.
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
      "-Ddraco_DIR=${draco_DIR}"
      "-DDRACO_LIBRARY_TYPE=${DRACO_LIBRARY_TYPE}"
      "-DDraco_TPL_INCLUDE_DIRS=${tmp}"
      "-DCMAKE_VERBOSE_MAKEFILE=TRUE"
      "-DDRACO_C4=${DRACO_C4}"
      )
    if( WIN32 )
      # For Win32 builds, DLL and applications are built in the directory
      # specified by CMAKE_RUNTIME_OUTPUT_DIRECTORY.
      set( build_system_state "${build_system_state};-DCMAKE_RUNTIME_OUTPUT_DIRECTORY=${Draco_BINARY_DIR}/\${CMAKE_BUILD_TYPE}" )
    else()
      set( build_system_state "${build_system_state};-DWITH_CUDA=${WITH_CUDA}" )
    endif()

    cmake_add_fortran_subdirectory(
      ftest                        # directory name
      PROJECT quadrature_ftest  # project name in toplevel CMakeLists.txt
      # Creates target named 'quadrature_ftest_build'
      ARCHIVE_DIR  ftest          # .lib location relative to root binary tree
      RUNTIME_DIR  ftest          # .dll location relative to root binary tree
      LIBRARIES    rtt_quadrature_ftest   # new library will have this name.
      # --> librtt_FC_f90.dll.
      TARGET_NAMES Lib_quadrature_ftest
      DEPENDS      "Lib_quadrature_test;Lib_quadrature;Lib_dsxx"
      NO_EXTERNAL_INSTALL
      CMAKE_COMMAND_LINE "${build_system_state}"
      # VERBOSE
      )

    # This directory contains c++ main code that wraps the Fortran tests.
    add_subdirectory( fctest )
  endif()

endif()

# ---------------------------------------------------------------------------- #
# Autodoc
# ---------------------------------------------------------------------------- #

process_autodoc_pages()

#------------------------------------------------------------------------------#
# End quadrature/CMakeLists.txt
#------------------------------------------------------------------------------#
