#-----------------------------*-cmake-*----------------------------------------#
# file   config/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for norms.
# note   Copyright (C) 2016, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.9.0)
project( norms CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB sources *.cc )
file( GLOB explicit_instantiations *_pt.cc )
file( GLOB headers *.hh )

# Make the header files available in the IDE.
if( MSVC_IDE OR ${CMAKE_GENERATOR} MATCHES Xcode )
   list( APPEND sources ${headers} )
endif()

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories(
   ${PROJECT_SOURCE_DIR}      # sources
   ${PROJECT_SOURCE_DIR}/..   # c4/ds++ header files
   ${dsxx_BINARY_DIR}         # ds++/config.h
   ${c4_BINARY_DIR}           # c4/config.h
   ${MPI_CXX_INCLUDE_PATH}
   )

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library(
   TARGET       Lib_norms
   TARGET_DEPS  Lib_c4
   LIBRARY_NAME ${PROJECT_NAME}
   SOURCES      "${sources}" )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_norms EXPORT draco-targets DESTINATION ${DBSCFGDIR}lib )
install( FILES ${headers} DESTINATION ${DBSCFGDIR}include/norms )

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
