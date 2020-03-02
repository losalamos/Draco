#-----------------------------*-cmake-*----------------------------------------#
# file   quadrature/ftest/CMakeLists.txt
# author Allan Wollaber
# date   Mon May 23 16:15:44 MDT 2016
# brief  Generate a Fortran testing library for quadrature
# note   Copyright (C) 2016-2019, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.9.0)
project( quadrature_ftest C Fortran )

# ---------------------------------------------------------------------------- #
# Setup - For Generators other than "Unix Makefiles", this folder is
# parsed independent of parent folders, so we need  some basic setup.
# ---------------------------------------------------------------------------- #
# Build system configuration files are located here.
if( NOT EXISTS ${draco_DIR} )
  message( FATAL_ERROR "can't find draco/config directory at ${draco_DIR}" )
else()
  set( CMAKE_MODULE_PATH ${draco_DIR} )
endif()

# Provide helper functions used by component CMakeLists.txt files.
include( component_macros )

# ---------------------------------------------------------------------------- #
# Build package libraries
# ---------------------------------------------------------------------------- #

# If we are using cmake_add_fortran_subdirectory to include this package, the
# library targets Lib_dsxx and Lib_quadrature will not be defined.  We must
# manually find these libraries.
#
# Note: A more complex example of CMakeAddFortranSubdirectory can be
# found in jayenne/src/api/ftest.

if(NOT TARGET Lib_dsxx)
  # Rebuild the list Draco_TPL_INCLUDE_DIRS from the packed list (see
  # api/CMakeLists.txt) by replacing triple underscores with a semicolon.
  # This must be done before calling find_package(draco)
  string( REGEX REPLACE "___" ";" cafs_Draco_TPL_INCLUDE_DIRS
    "${Draco_TPL_INCLUDE_DIRS}" )

  include(CMakeAddFortranSubdirectory)

  get_filename_component( draco_BINARY_DIR ${PROJECT_BINARY_DIR}/../../.. 
    ABSOLUTE )
  cafs_create_imported_targets( Lib_dsxx
    "rtt_ds++"       "${draco_BINARY_DIR}/src/ds++"       CXX )
  cafs_create_imported_targets( Lib_c4 "rtt_c4"
    "${draco_BINARY_DIR}/src/c4" CXX )
  cafs_create_imported_targets( Lib_quadrature
    "rtt_quadrature" "${draco_BINARY_DIR}/src/quadrature" CXX )
  cafs_create_imported_targets( Lib_quadrature_test
    "rtt_quadrature_test" "${draco_BINARY_DIR}/src/quadrature/test" CXX )

  # If we get here, we also need to use the Draco scripts to setup compiler
  # flags and MPI options
  include( buildEnv )
  dbsSetDefaults()
  # On Win32, set the default top level output directory:
  if(WIN32)
    set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${draco_BINARY_DIR}/${CMAKE_BUILD_TYPE}
      CACHE PATH "Build runtime objects at this location.")
  endif()
  include( platform_checks )
  include( compilerEnv )
  dbsSetupFortran()
  # GSL is needed by quadrature.
  include( vendor_libraries )
  setupVendorLibraries()

  # link to msys64 formatted msmpi.a instead of VS formatted msmpi.lib/dll,
  # and ensure that the compiler option '-fno-range-check' is enabled.
  cafs_fix_mpi_library()

endif()

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #
set( f90sources
  ${PROJECT_SOURCE_DIR}/tstquadrature_interfaces.f90
  ${PROJECT_SOURCE_DIR}/../quadrature_interfaces.f90 )

# ---------------------------------------------------------------------------- #
# Build library for test directory
# ---------------------------------------------------------------------------- #

# Xcode: Since cafs_create_imported_targets does not set all of the 
# dependencies between the imported libraries, we need to list them explicitly.
set(target_deps
    Lib_quadrature_test
    Lib_quadrature
    Lib_dsxx )
if( "${DRACO_C4}" STREQUAL "MPI" )
  add_definitions( -DC4_MPI )
  list(APPEND target_deps MPI::MPI_Fortran )
endif()

add_component_library(
  TARGET       Lib_quadrature_ftest
  LIBRARY_NAME quadrature_ftest
  SOURCES      "${f90sources}"
  LINK_LANGUAGE "Fortran"
  TARGET_DEPS   "${target_deps}"
  NOEXPORT )

#------------------------------------------------------------------------------#
# end quadrature/ftest/CMakeLists.txt
#------------------------------------------------------------------------------#
