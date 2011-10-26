#-----------------------------*-cmake-*----------------------------------------#
# file   device/GPU_Device.cmake
# author Kelly Thompson
# date   2011 Oct 24
# brief  Instructions for building device Makefile.
# note   © Copyright 2011 Los Alamos National Security, All rights reserved.
#------------------------------------------------------------------------------#
# $Id$

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #

configure_file( config.h.in ${PROJECT_BINARY_DIR}/device/config.h )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( sources 
   GPU_Device.cc
   )
set( headers 
   GPU_Device.hh
   ${PROJECT_BINARY_DIR}/device/config.h
   )
# file( GLOB cudaSources *.cu )

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
   ${CUDA_TOOLKIT_INCLUDE}
   )

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library( Lib_device device "${sources}" )
target_link_libraries( Lib_device 
   Lib_dsxx
   ${CUDA_CUDA_LIBRARY} 
   )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_device DESTINATION lib )
install( FILES ${headers} DESTINATION include/device )

