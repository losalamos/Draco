#--------------------------------------------*-cmake-*---------------------------------------------#
# file   c4/ftest/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2013 June 24
# brief  Generate build project files for Fortran.
# note   Copyright (C) 2016-2020, Triad National Security, LLC.
#        All rights reserved.
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.17.0)
project( c4_ftest C Fortran )

#--------------------------------------------------------------------------------------------------#
# If this directory is built via CAFS (ie. MSVC or XCode), then we need to do
# some extra setup.  When build via CAFS, target Lib_dsxx won't be defined yet.
#--------------------------------------------------------------------------------------------------#
if(NOT TARGET Lib_dsxx)
  include(CMakeAddFortranSubdirectory)
  CAFS_Fortran_dir_boilerplate_setup()
  CAFS_Fortran_dir_MPI_setup()

  # CAFS setup unique to this directory.
  # 1. create import targets that we will link against.
  if( NOT fc4_BINARY_DIR )
    set( fc4_BINARY_DIR "${draco_BINARY_DIR}/src/c4/fc4" )
  endif()
  cafs_create_imported_targets(Lib_c4_fc4 "librtt_fc4" "${fc4_BINARY_DIR}" CXX)
endif()

#--------------------------------------------------------------------------------------------------#
# Source files
#--------------------------------------------------------------------------------------------------#
set( f90sources
  ${PROJECT_SOURCE_DIR}/fc4_test.f90
  ${PROJECT_SOURCE_DIR}/mpi_hw_ftest.f90 )

#--------------------------------------------------------------------------------------------------#
# Build library for test directory
#--------------------------------------------------------------------------------------------------#
set(target_deps Lib_c4_fc4 Lib_c4 Lib_dsxx )
set(include_dirs PUBLIC
  "$<BUILD_INTERFACE:${c4_ftest_BINARY_DIR}>"
  "$<BUILD_INTERFACE:${fc4_BINARY_DIR}>")
add_component_library(
  TARGET        Lib_c4_ftest
  LIBRARY_NAME  c4_ftest
  SOURCES       "${f90sources}"
  LINK_LANGUAGE Fortran
  TARGET_DEPS   "${target_deps}"
  INCLUDE_DIRS  "${include_dirs}"
  NOEXPORT )

#--------------------------------------------------------------------------------------------------#
# End c4/ftest/CMakeLists.txt
#--------------------------------------------------------------------------------------------------#
