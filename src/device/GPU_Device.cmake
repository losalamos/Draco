#-----------------------------*-cmake-*----------------------------------------#
# file   device/GPU_Device.cmake
# brief  Instructions for building device Makefile.
# note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$

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
   GPU_Module.cc
   )
set( headers 
   GPU_Device.hh
   GPU_Module.hh
   ${PROJECT_BINARY_DIR}/device/config.h
   device_cuda.h
   )

# Make the header files available in the IDE.
if( MSVC_IDE OR ${CMAKE_GENERATOR} MATCHES Xcode)
   list( APPEND sources ${headers} )
endif()

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories(
   ${PROJECT_SOURCE_DIR}        # sources
   ${PROJECT_BINARY_DIR}        # config.h
   ${draco_src_dir_SOURCE_DIR}  # ds++ header files
   ${dsxx_BINARY_DIR}           # ds++/config.h
   )

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library( 
   TARGET       Lib_device
   TARGET_DEPS  Lib_dsxx
   LIBRARY_NAME device 
   SOURCES      "${sources}"
   VENDOR_LIST  "CUDA"
   VENDOR_LIBS  "${CUDA_CUDA_LIBRARY}"
   VENDOR_INCLUDE_DIRS "${CUDA_TOOLKIT_INCLUDE}"
   )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_device EXPORT draco-targets DESTINATION ${DBSCFGDIR}lib )
install( FILES ${headers} DESTINATION ${DBSCFGDIR}include/device )
