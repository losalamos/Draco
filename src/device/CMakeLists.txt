#-----------------------------*-cmake-*----------------------------------------#
# file   device/CMakeLists.txt
# author Gabriel Rockefeller
# date   2011 June 13
# brief  Instructions for building device Makefile.
# note   © Copyright 2011 Los Alamos National Security, All rights reserved.
#------------------------------------------------------------------------------#
# $Id$

cmake_minimum_required( VERSION 2.8 FATAL_ERROR )
project( device CXX )

if( NOT EXISTS ${CUDA_NVCC_EXECUTABLE} )
   find_package( CUDA QUIET )
   # find_package will return CUDA_FOUND=TRUE if libcuda is
   # available.  We also need the nvcc executable
endif()

# ---------------------------------------------------------------------------- #
# Build one of DACS_Device or GPU_Device
# ---------------------------------------------------------------------------- #

if( ENABLE_HET_CODE )
   include( DACS_Device.cmake )
elseif( USE_CUDA AND EXISTS "${CUDA_NVCC_EXECUTABLE}" )
   include( GPU_Device.cmake )
endif()

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

# Push some variables up one level
set( Draco_EXPORT_TARGET_PROPERTIES
   "${Draco_EXPORT_TARGET_PROPERTIES}" PARENT_SCOPE)
