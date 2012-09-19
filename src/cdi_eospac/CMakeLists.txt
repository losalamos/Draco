#-----------------------------*-cmake-*----------------------------------------#
# file   src/cdi_eospac/CMakeLists.txt
# brief  Instructions for building cdi_eospac Makefiles.
# note   Copyright (C) 2010-2012 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id: CMakeLists.txt 6732 2012-09-05 22:28:18Z kellyt $
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 2.6)
project( cdi_eospac CXX )

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #

# Create cdi_eospac/config.h
set(V_INT   int)
set(V_FLOAT double)
set(V_BOOL  int)
configure_file( config.h.in ${PROJECT_BINARY_DIR}/cdi_eospac/config.h )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB sources *.cc )
file( GLOB headers *.hh )
list( APPEND headers ${PROJECT_BINARY_DIR}/cdi_eospac/config.h )

# Make the header files available in the IDE.
if( MSVC_IDE OR ${CMAKE_GENERATOR} MATCHES Xcode )
   list( APPEND sources ${headers} )
endif()

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories( 
   ${PROJECT_SOURCE_DIR}      # component headers
   ${PROJECT_BINARY_DIR}      # cdi_eospac/config.h
   ${PROJECT_SOURCE_DIR}/..   # ds++ and cdi header files
   ${dsxx_BINARY_DIR}         # ds++/config.h
   ${ode_BINARY_DIR}          # ode/config.h
#    ${parser_BINARY_DIR}
) 

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library( 
   TARGET       Lib_cdi_eospac
   TARGET_DEPS  Lib_cdi
   LIBRARY_NAME cdi_eospac
   SOURCES      "${sources}"  )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_cdi_eospac EXPORT draco-targets DESTINATION lib )
install( FILES ${headers} DESTINATION include/cdi_eospac )

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
