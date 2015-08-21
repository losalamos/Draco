#-----------------------------*-cmake-*----------------------------------------#
# file   FortranChecks/f90sub/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2013 June 24
# brief  Generate build project files for fortran.
# note   Copyright (C) 2013-2015, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id: CMakeLists.txt 6705 2012-08-29 16:21:06Z kellyt $
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.0.0)
project( FortranChecks_f90sub Fortran )

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
# library targets Lib_dsxx and Lib_FC_derived_typ will not be defined.
# We must manually find these libraries.
#
# Note: A more complex example of CMakeAddFortranSubdirectory can be
# found in jayenne/src/wedgehog/ftest.

if(NOT TARGET Lib_dsxx)
  # Rebuild the list Draco_TPL_INCLUDE_DIRS from the packed list (see wedgehog/CMakeLists.txt)
  # by replacing triple underscores with a semicolon.  This must be done before calling
  # find_package(draco)
  string( REGEX REPLACE "___" ";" cafs_Draco_TPL_INCLUDE_DIRS "${Draco_TPL_INCLUDE_DIRS}" )

  include(CMakeAddFortranSubdirectory)
  get_filename_component( draco_BINARY_DIR ${PROJECT_BINARY_DIR}/../../.. ABSOLUTE )
  cafs_create_imported_targets( Lib_dsxx "rtt_ds++" "${draco_BINARY_DIR}/src/ds++" CXX )
  cafs_create_imported_targets( Lib_FC_Derived_Type "rtt_FC_Derived_Type" "${draco_BINARY_DIR}/src/FortranChecks" CXX )

  # If we get here, we also need to use the Draco scripts to setup compiler flags and MPI options
  include( buildEnv )
  dbsSetDefaults()
  # On Win32, set the default top level output directory:
  if(WIN32)
    set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${draco_BINARY_DIR}/${CMAKE_BUILD_TYPE}
      CACHE PATH "Build runtime objects at this location.")
  endif()
  include( platform_checks )
  query_openmp_availability()
  include( compilerEnv )
  dbsSetupFortran()
  #include( vendor_libraries )
  #setupVendorLibraries()
endif()

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #
file( GLOB f90sources *.f90 )

add_component_library(
  TARGET       Lib_FC_f90sub
  LIBRARY_NAME FC_f90sub
  SOURCES      "${f90sources}"
  LINK_LANGUAGE "Fortran"
  TARGET_DEPS   "Lib_FC_Derived_Type;Lib_dsxx"
  NOEXPORT )

# Debug target properties:
#include(print_target_properties)
#echo_targets("Lib_FC_f90sub;Lib_FC_Derived_Type")
