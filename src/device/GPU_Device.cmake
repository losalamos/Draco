#-----------------------------*-cmake-*----------------------------------------#
# file   device/GPU_Device.cmake
# brief  Instructions for building device Makefile.
# note   Copyright (C) 2016-2019 Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

# Ref: https://devblogs.nvidia.com/building-cuda-applications-cmake/
#      https://github.com/robertmaynard/code-samples/blob/master/posts/cmake
# Ref: "Acceleware CUDA Course Lectures.pdf"

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #

set( TEST_KERNEL_BINDIR ${PROJECT_BINARY_DIR}/test CACHE PATH
  "GPU kernel binary install location" )
configure_file( config.h.in ${PROJECT_BINARY_DIR}/device/config.h )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( sources
   GPU_Device.cc
   GPU_Module.cc )
set( headers
   GPU_Device.hh
   GPU_Module.hh
   ${PROJECT_BINARY_DIR}/device/config.h
   device_cuda.h )

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library(
   TARGET       Lib_device
   TARGET_DEPS  Lib_dsxx
   LIBRARY_NAME device
   SOURCES      "${sources}"
   HEADERS      "${headers}" )
target_include_directories( Lib_device
  PUBLIC $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}> )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_device EXPORT draco-targets DESTINATION ${DBSCFGDIR}lib )
install( FILES ${headers} DESTINATION ${DBSCFGDIR}include/device )

#------------------------------------------------------------------------------#
# End device/GPU_Device.cmake
#------------------------------------------------------------------------------#
