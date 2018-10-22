#-----------------------------*-cmake-*----------------------------------------#
# file   FortranChecks/f90sub/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2013 June 24
# brief  Generate build project files for fortran.
# note   Copyright (C) 2016-2018, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.9.0)
project( FortranChecks_f90sub Fortran )

# ---------------------------------------------------------------------------- #
# Setup - For Generators other than "Unix Makefiles", this folder is parsed
# independent of parent folders, so we need some basic setup.
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
# library targets Lib_dsxx and Lib_FC_derived_typ will not be defined.  We must
# manually find these libraries.
#
# Note: A more complex example of CMakeAddFortranSubdirectory can be found in
# jayenne/src/wedgehog/ftest.

if(NOT TARGET Lib_dsxx)
  # Rebuild the list Draco_TPL_INCLUDE_DIRS from the packed list (see
  # wedgehog/CMakeLists.txt) by replacing triple underscores with a semicolon.
  # This must be done before calling find_package(draco)
  string( REGEX REPLACE "___" ";" cafs_Draco_TPL_INCLUDE_DIRS
    "${Draco_TPL_INCLUDE_DIRS}" )

  include(CMakeAddFortranSubdirectory)

  get_filename_component( draco_BINARY_DIR ${PROJECT_BINARY_DIR}/../../..
    ABSOLUTE )
  cafs_create_imported_targets( Lib_dsxx "rtt_ds++"
    "${draco_BINARY_DIR}/src/ds++" CXX )
  cafs_create_imported_targets( Lib_FC_Derived_Type "rtt_FC_Derived_Type"
    "${draco_BINARY_DIR}/src/FortranChecks" CXX )

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
  # Vendors are not needed here, but I'm trying to follow the model found in
  # wedgehog/ftest.  I've included these here to ensure that the slave build can
  # parse the find vendor scripts.
  include( vendor_libraries )
  setupVendorLibraries()

  # link to msys64 formatted msmpi.a instead of VS formatted msmpi.lib/dll,
  # and ensure that the compiler option '-fno-range-check' is enabled.
  cafs_fix_mpi_library()

endif()

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #
file( GLOB f90sources *.f90 *.F90 )

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #
if( DEFINED cafs_Draco_TPL_INCLUDE_DIRS )
  include_directories( "${cafs_Draco_TPL_INCLUDE_DIRS}" )
elseif( DEFINED MPI_Fortran_INCLUDE_PATH )
  # Only include directories if mpif.h is found.
  set( mpifh_found FALSE )
  foreach( dir ${MPI_Fortran_INCLUDE_PATH} )
    if( EXISTS ${dir}/mpif.h )
      set( mpifh_found TRUE )
    endif()
  endforeach()
  if( mpifh_found )
    include_directories( "${MPI_Fortran_INCLUDE_PATH}" )
  endif()
endif()
if( "${DRACO_C4}" STREQUAL "MPI" )
  add_definitions( -DC4_MPI )
endif()

# ---------------------------------------------------------------------------- #
# Build library for test directory
# ---------------------------------------------------------------------------- #

# Xcode: Since cafs_create_imported_targets does not set all of the dependencies
# between the imported libraries, we need to list them explicitly.
set(target_deps "Lib_FC_Derived_Type;Lib_dsxx;${MPI_Fortran_LIBRARIES}" )

add_component_library(
  TARGET       Lib_FC_f90sub
  LIBRARY_NAME FC_f90sub
  SOURCES      "${f90sources}"
  LINK_LANGUAGE "Fortran"
  TARGET_DEPS   "${target_deps}"
  NOEXPORT )

#------------------------------------------------------------------------------#
# End f90sub/CMakeLists.txt
#------------------------------------------------------------------------------#
