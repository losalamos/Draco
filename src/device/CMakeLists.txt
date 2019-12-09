#-----------------------------*-cmake-*----------------------------------------#
# file   device/CMakeLists.txt
# author Gabriel Rockefeller
# date   2011 June 13
# brief  Instructions for building device Makefile.
# note   Copyright (C) 2016-2019 Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.11.0)
project( device CXX ${CUDA_DBS_STRING})

# ---------------------------------------------------------------------------- #
# Build one of DACS_Device or GPU_Device
# ---------------------------------------------------------------------------- #

if( WITH_CUDA )
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

# ---------------------------------------------------------------------------- #
# End devices/CMakeLists.txt
# ---------------------------------------------------------------------------- #
