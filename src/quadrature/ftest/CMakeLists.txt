#-----------------------------*-cmake-*----------------------------------------#
# file   quadrature/ftest/CMakeLists.txt
# author Allan Wollaber
# date   Mon May 23 16:15:44 MDT 2016
# brief  Generate a Fortran testing library for quadrature
# note   Copyright (C) 2016, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.9.0)
project( quadrature_ftest Fortran )

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
# found in jayenne/src/wedgehog/ftest.

if(NOT TARGET Lib_dsxx)
  # Rebuild the list Draco_TPL_INCLUDE_DIRS from the packed list (see
  # wedgehog/CMakeLists.txt) by replacing triple underscores with a semicolon.
  # This must be done before calling find_package(draco)
  string( REGEX REPLACE "___" ";" cafs_Draco_TPL_INCLUDE_DIRS
    "${Draco_TPL_INCLUDE_DIRS}" )

  # message("dtplid = ${cafs_Draco_TPL_INCLUDE_DIRS}")

  include(CMakeAddFortranSubdirectory)
  get_filename_component( draco_BINARY_DIR ${PROJECT_BINARY_DIR}/../../.. ABSOLUTE )
  cafs_create_imported_targets( Lib_dsxx       
    "rtt_ds++"       "${draco_BINARY_DIR}/src/ds++"       CXX )
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
  #query_openmp_availability()
  include( compilerEnv )
  dbsSetupFortran()
  # GSL is needed by quadrature.
  include( vendor_libraries )
  setupVendorLibraries()

  if(WIN32)
    # MS-MPI 5 and gfortran/4.8 do not play nice together...
    if( "${CMAKE_BUILD_TYPE}x" STREQUAL "x")
      message( FATAL_ERROR
        "In ftest/CMakeLists.txt, CMAKE_BUILD_TYPE must be set."
        "Note: for multi-config build tools, CMakeAddFortranSubdirectory will "
        "set this for you.")
    endif()
    string( TOUPPER ${CMAKE_BUILD_TYPE} ucbt )
    # TODO: Can this be removed?
    if( "${CMAKE_Fortran_FLAGS_${ucbt}}" MATCHES "frange-check" )
      string(REPLACE "range-check" "no-range-check" CMAKE_Fortran_FLAGS_DEBUG
        ${CMAKE_Fortran_FLAGS_DEBUG} )
    else()
      set( CMAKE_Fortran_FLAGS_${ucbt}
        "${CMAKE_Fortran_FLAGS_${ucbt}} -fno-range-check")
    endif()

    # Don't link to the C++ MS-MPI library when compiling with MinGW
    # gfortran. Instead, link to libmsmpi.a that was created via gendef.exe and
    # dlltool.exe from msmpi.dll.  Ref:
    # http://www.geuz.org/pipermail/getdp/2012/001519.html
    if( "${MPI_Fortran_LIBRARIES}none" STREQUAL "none" OR
        "${MPI_Fortran_LIBRARIES}" MATCHES "msmpi.lib" )
      # should be located in ENV{PATH} at c:/Program Files/Microsoft MPI/Bin/
      find_file( MPI_gfortran_LIBRARIES NAMES "libmsmpi.a" )
      set( MPI_Fortran_LIBRARIES ${MPI_gfortran_LIBRARIES} CACHE FILEPATH
        "msmpi for gfortran" FORCE )
    endif()
  endif()

endif()

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #
set( f90sources
  ${PROJECT_SOURCE_DIR}/tstquadrature_interfaces.f90
  ${PROJECT_SOURCE_DIR}/../quadrature_interfaces.f90
   )

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
set(target_deps 
    Lib_quadrature_test
    Lib_quadrature
    Lib_dsxx
    "${MPI_Fortran_LIBRARIES}" )

add_component_library(
  TARGET       Lib_quadrature_ftest
  LIBRARY_NAME quadrature_ftest
  SOURCES      "${f90sources}"
  LINK_LANGUAGE "Fortran"
  TARGET_DEPS   "${target_deps}"
  NOEXPORT )

# Debug target properties:
#include(print_target_properties)
#print_targets_properties( "Lib_dsxx;Lib_quadrature;Lib_quadrature_test" )

#------------------------------------------------------------------------------#
# end quadrature/ftest/CMakeLists.txt
#------------------------------------------------------------------------------#
