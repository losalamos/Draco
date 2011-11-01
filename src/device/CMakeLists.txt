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
endif()

# ---------------------------------------------------------------------------- #
# Build one of DACS_Device or GPU_Device
# ---------------------------------------------------------------------------- #

if( ENABLE_HET_CODE )
   include( DACS_Device.cmake )
elseif( EXISTS "${CUDA_CUDA_LIBRARY}" )
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


