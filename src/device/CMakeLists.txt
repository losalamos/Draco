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

# At the moment, this package is het-only, and x86-only, but it builds
# a small ppe binary for testing (so that dacs_de_start can have
# something to launch).

if( NOT "${CMAKE_CXX_COMPILER}" MATCHES "[sp]pu-g[+][+]" )  # if x86 build

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #

   if( NOT EXISTS ${PPE_PREFIX}/bin/dacs_device_ppe_exe )
      string ( REPLACE "x86" "ppe" PPE_PREFIX ${CMAKE_INSTALL_PREFIX} )
      set( PPE_PREFIX ${PPE_PREFIX} CACHE PATH 
         "PPE binary install location" )
   endif()

   # sanity check
   if ( NOT EXISTS ${PPE_PREFIX}/bin/dacs_device_ppe_exe )
      message(FATAL_ERROR "PPE_PREFIX must be set to the installed PPE"
         " build of Draco." )
   endif()

   include(draco_macros)
   draco_flat_function_names()
   configure_file( config.h.in ${PROJECT_BINARY_DIR}/device/config.h )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

   file( GLOB template_implementations *.t.hh *.i.hh )
   file( GLOB sources *.cc )
   file( GLOB explicit_instantiations *_pt.cc )
   file( GLOB headers *.hh )
   list( APPEND headers ${PROJECT_BINARY_DIR}/device/config.h )

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

   add_component_library( Lib_device device "${sources}" )
   target_link_libraries( Lib_device 
      /usr/lib64/libdacs_hybrid.so
      )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

   install( TARGETS Lib_device DESTINATION lib )
   install( FILES ${headers} DESTINATION include/device )

endif() # endif x86 build

# ---------------------------------------------------------------------------- #
# Unit tests
# ---------------------------------------------------------------------------- #

if( BUILD_TESTING )
   add_subdirectory( test )
endif()   
