#-----------------------------*-cmake-*----------------------------------------#
# file   config/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for fpe_trap.
# note   Copyright (C) 2010-2012, Los Alamos National Security
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 2.6)
project( fpe_trap CXX )

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #

if( ${CMAKE_SYSTEM} MATCHES "OSF" )
   set( FPETRAP_OSF_ALPHA 1 )
elseif( ${CMAKE_SYSTEM} MATCHES "Linux" )
   if( CMAKE_COMPILER_IS_GNUCXX )
      set( FPETRAP_LINUX_X86 1 )
   endif()
elseif( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
   if( ${_CMAKE_OSX_MACHINE} MATCHES "ppc" )
      set( FPETRAP_DARWIN_PPC 1 )
   else()
      set( FPETRAP_DARWIN_INTEL 1 )
   endif()
endif()

configure_file( config.h.in ${PROJECT_BINARY_DIR}/fpe_trap/config.h )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB sources *.cc )
file( GLOB explicit_instantiations *_pt.cc )
file( GLOB headers *.hh )
list( APPEND headers ${PROJECT_BINARY_DIR}/fpe_trap/config.h )

# Make the header files available in the IDE.
if( MSVC_IDE OR ${CMAKE_GENERATOR} MATCHES Xcode )
   list( APPEND sources ${headers} )
endif()

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories( 
   ${PROJECT_SOURCE_DIR}      # sources
   ${PROJECT_BINARY_DIR}      # config.h
   ${PROJECT_SOURCE_DIR}/..   # ds++ header files
   ${dsxx_BINARY_DIR}         # ds++/config.h
) 

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library( 
   TARGET       Lib_fpe_trap 
   TARGET_DEPS  Lib_dsxx
   LIBRARY_NAME ${PROJECT_NAME} 
   SOURCES      "${sources}"  )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_fpe_trap EXPORT draco-targets DESTINATION lib )
install( FILES ${headers} DESTINATION include/fpe_trap )

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


