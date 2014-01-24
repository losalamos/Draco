#-----------------------------*-cmake-*----------------------------------------#
# file   device/DACS_Device.cmake
# author Gabriel Rockefeller
# date   2011 June 13
# brief  Instructions for building device Makefile.
# note   Copyright (C) 2011-2014 Los Alamos National Security, 
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

# At the moment, this package is het-only, and x86-only, but it builds
# small ppe binaries for testing (so that dacs_de_start can have
# something to launch).

# BUG in <dacs.h>: pgCC (with strict ansi flag -A) will not compile
# code that includes <dacs.h> because of trailing commas found in
# enumeration lists.

if( NOT "${CMAKE_CXX_COMPILER}" MATCHES "[sp]pu-g[+][+]" )  # if x86 build

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #

   add_feature_info( TEST_PPE_BINDIR TEST_PPE_BINDIR
      "Provide full path to the cell binary dacs_noop_ppe_exe when it
      is not in the default location (x86 and ppe directories in
      parallel locations). ")

   # Record the location of the Cell-side test binaries in config.h.
   if( NOT EXISTS ${TEST_PPE_BINDIR}/dacs_noop_ppe_exe )
      string ( REPLACE "x86" "ppe" TEST_PPE_BINDIR 
         "${CMAKE_INSTALL_PREFIX}/bin" )
      # Simply path name
      get_filename_component( TEST_PPE_BINDIR ${TEST_PPE_BINDIR}
         ABSOLUTE )
      # Save it in the cache file.
      set( TEST_PPE_BINDIR ${TEST_PPE_BINDIR} CACHE PATH 
         "PPE binary install location" )
   endif()

   # Sanity check
   if ( NOT EXISTS ${TEST_PPE_BINDIR}/dacs_noop_ppe_exe )
      message(FATAL_ERROR "TEST_PPE_BINDIR must be set to bin "
         "directory of the installed PPE build of Draco." )
   endif()
   
   # Simply path name
   get_filename_component( TEST_PPE_BINDIR ${TEST_PPE_BINDIR}
      ABSOLUTE )
   configure_file( config.h.in ${PROJECT_BINARY_DIR}/device/config.h )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

   set( sources 
      DACS_Device.cc
      DACS_Device_Interface.cc
      DACS_Process.cc
      )
   set( headers 
      DACS_Device.hh
      DACS_Device_Interface.hh
      DACS_External_Process.hh
      DACS_Process.hh
      ${PROJECT_BINARY_DIR}/device/config.h
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
      VENDOR_LIST  "DaCS"
      VENDOR_LIBS  "/opt/ofed/lib64/librdmacm.so;/opt/ofed/lib64/libibverbs.so;/opt/PBS/lib64/libtorque.so;/usr/lib64/libnuma.so;/usr/lib64/libdacs_hybrid.so"
      )
# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

   install( TARGETS Lib_device EXPORT draco-targets DESTINATION lib )
   install( FILES ${headers} DESTINATION include/device )

endif() # endif x86 build
