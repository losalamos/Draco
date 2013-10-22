#-----------------------------*-cmake-*----------------------------------------#
# file   ds++/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 Sep 24
# brief  Generate build project files for ds++.
# note   Copyright (C) 2010-2012, Los Alamos National Security
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 2.6)
project( dsxx CXX )

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #

# Review config/platform_checks.cmake for many of these settings.

# Choose platforms for fpe_trap
if( ${CMAKE_SYSTEM} MATCHES "OSF" )
   set( FPETRAP_OSF_ALPHA 1 )
elseif( ${CMAKE_SYSTEM} MATCHES "Linux" )
   # If the OS is Linux, fpe_trap should work as long as
   # feenableexscept is available.  That is, it should work for both
   # GNU and Intel compiler suites.
   #if( CMAKE_COMPILER_IS_GNUCXX )
      set( FPETRAP_LINUX_X86 1 )
   #endif()
   # Turn fpe_trap off for PGI compilers.
   if( ${CMAKE_C_COMPILER} MATCHES pgcc )
      set( FPETRAP_LINUX_X86 0 )
   endif()
elseif( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
   if( ${_CMAKE_OSX_MACHINE} MATCHES "ppc" )
      set( FPETRAP_DARWIN_PPC 1 )
   else()
      set( FPETRAP_DARWIN_INTEL 1 )
   endif()
elseif( ${CMAKE_SYSTEM_NAME} MATCHES "Windows" )
   set( FPETRAP_WINDOWS_X86 1 )
endif()

# Create ds++/config.h
configure_file( config.h.in ${PROJECT_BINARY_DIR}/ds++/config.h )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB template_implementations *.t.hh *.i.hh )
file( GLOB sources *.cc )
file( GLOB explicit_instantiations *_pt.cc )
file( GLOB headers *.hh )
list( APPEND headers ${PROJECT_BINARY_DIR}/ds++/config.h )
list( REMOVE_ITEM headers ${template_implementations} )

# BUG: artf21615 (tf.lanl.gov): Safe_Ptr has a design defect
list( REMOVE_ITEM headers 
   ${PROJECT_SOURCE_DIR}/Safe_Ptr.hh
   ${PROJECT_SOURCE_DIR}/DBC_Ptr.hh 
)

# Make the header files available in the IDE.
if( MSVC_IDE OR ${CMAKE_GENERATOR} MATCHES Xcode )
   list( APPEND sources ${headers} ${template_implementations} )
endif()

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories( 
   ${PROJECT_SOURCE_DIR}   # sources
   ${PROJECT_BINARY_DIR} ) # config.h

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library( 
   TARGET       Lib_dsxx 
   LIBRARY_NAME ds++ 
   SOURCES      "${sources}" 
   )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_dsxx EXPORT draco-targets DESTINATION lib )
target_link_libraries( Lib_dsxx ${Lib_win_winsock} )
install( FILES ${headers} ${template_implementations} DESTINATION include/ds++ )

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
