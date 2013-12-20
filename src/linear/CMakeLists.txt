#-----------------------------*-cmake-*----------------------------------------#
# file   config/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for linear.
# note   Copyright (C) 2010-2012, Los Alamos National Security
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#
cmake_minimum_required( VERSION 2.8 )
project( linear CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB template_implementations *.t.hh *.i.hh )
file( GLOB sources *.cc )
file( GLOB headers *.hh )
list( REMOVE_ITEM headers ${template_implementations} )

# Make the header files available in the IDE.
if( MSVC_IDE OR ${CMAKE_GENERATOR} MATCHES Xcode )
   list( APPEND sources ${headers} )
endif()

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories( 
   ${PROJECT_SOURCE_DIR}      # sources
   ${PROJECT_BINARY_DIR}      # linear/config.h
   ${PROJECT_SOURCE_DIR}/..   # c4/ds++ header files
   ${dsxx_BINARY_DIR}         # ds++/config.h
) 

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library(
   TARGET       Lib_linear
   TARGET_DEPS  Lib_c4
   LIBRARY_NAME ${PROJECT_NAME}
   SOURCES      "${sources}" )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_linear EXPORT draco-targets DESTINATION lib )
install( FILES ${headers} ${template_implementations} DESTINATION include/linear )

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
