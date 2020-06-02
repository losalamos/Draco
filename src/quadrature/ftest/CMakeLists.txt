#-----------------------------*-cmake-*----------------------------------------#
# file   quadrature/ftest/CMakeLists.txt
# author Allan Wollaber
# date   Mon May 23 16:15:44 MDT 2016
# brief  Generate a Fortran testing library for quadrature
# note   Copyright (C) 2016-2020, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.17.0)
project( quadrature_ftest C Fortran )

#------------------------------------------------------------------------------#
# If this directory is built via CAFS (ie. MSVC or XCode), then we need to do
# some extra setup.  When build via CAFS, target Lib_dsxx won't be defined yet.
#------------------------------------------------------------------------------#
if(NOT TARGET Lib_dsxx)
  include(CMakeAddFortranSubdirectory)
  CAFS_Fortran_dir_boilerplate_setup()
  CAFS_Fortran_dir_MPI_setup() # returns CAFS_MPI_DEPS

  # CAFS setup unique to this directory.
  # 1. create import targets that we will link against.
  cafs_create_imported_targets( Lib_quadrature "rtt_quadrature"
    "${draco_BINARY_DIR}/src/quadrature" CXX )
  cafs_create_imported_targets( Lib_quadrature_test "rtt_quadrature_test"
    "${draco_BINARY_DIR}/src/quadrature/test" CXX )
  # 2. Vendor discovery is not needed for this package
  # setupVendorLibraries()
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
set(target_deps Lib_quadrature_test Lib_quadrature Lib_dsxx ${CAFS_MPI_DEPS} )
#if( "${DRACO_C4}" STREQUAL "MPI" )
#  list(APPEND target_deps MPI::MPI_Fortran )
#endif()
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
