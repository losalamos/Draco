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
# small ppe binaries for testing (so that dacs_de_start can have
# something to launch).

# BUG in <dacs.h>: pgCC (with strict ansi flag -A) will not compile
# code that includes <dacs.h> because of trailing commas found in
# enumeration lists.  I have provided a 'fixed' version of dacs.h in
# this directory that should be used in place of the system installed
# version found at /opt/cell/sysroot/usr/include/.

if( NOT "${CMAKE_CXX_COMPILER}" MATCHES "[sp]pu-g[+][+]" )  # if x86 build

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #

   # Record the location of the Cell-side test binaries in config.h.
   if( NOT EXISTS ${TEST_PPE_BINDIR}/dacs_noop_ppe_exe )
      string ( REPLACE "x86" "ppe" TEST_PPE_BINDIR 
         "${CMAKE_INSTALL_PREFIX}/bin" )
      set( TEST_PPE_BINDIR ${TEST_PPE_BINDIR} CACHE PATH 
         "PPE binary install location" )
   endif()

   # Sanity check
   if ( NOT EXISTS ${TEST_PPE_BINDIR}/dacs_noop_ppe_exe )
      message(FATAL_ERROR "TEST_PPE_BINDIR must be set to bin "
         "directory of the installed PPE build of Draco." )
   endif()

   configure_file( config.h.in ${PROJECT_BINARY_DIR}/device/config.h )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

   file( GLOB template_implementations *.t.hh *.i.hh )
   file( GLOB sources *.cc )
   file( GLOB explicit_instantiations *_pt.cc )
   file( GLOB headers *.hh *.h )
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
  

# ---------------------------------------------------------------------------- #
# Autodoc
# ---------------------------------------------------------------------------- #

process_autodoc_pages()


