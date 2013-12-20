#-----------------------------*-cmake-*----------------------------------------#
# file   config/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for memory.
# note   Copyright (C) 2010-2013, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id: CMakeLists.txt 7239 2013-10-07 20:29:39Z kellyt $
#------------------------------------------------------------------------------#
cmake_minimum_required( VERSION 2.8 )
project( memory CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB template_implementations *.t.hh *.i.hh )
file( GLOB sources *.cc )
file( GLOB explicit_instantiations *_pt.cc )
file( GLOB headers *.hh )
list( REMOVE_ITEM headers
   ${template_implementations} 
   ${PROJECT_SOURCE_DIR}/draco_info.cc
   )

# Make the header files available in the IDE.
if( MSVC_IDE OR ${CMAKE_GENERATOR} MATCHES Xcode )
   list( APPEND sources ${headers} )
endif()

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories(
   ${PROJECT_SOURCE_DIR}  # sources
   ${PROJECT_BINARY_DIR}  # config.h
   ${draco_src_dir_SOURCE_DIR} # ds++ and c4 headers
   ${dsxx_BINARY_DIR} # ds++/config.h
   ${diagnostics_BINARY_DIR} $ diagnostics/config.h
   ${c4_BINARY_DIR}   # c4/config.h
   ${MPI_CXX_INCLUDE_PATH}
   )

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library( 
   TARGET       Lib_memory
   TARGET_DEPS  Lib_c4
   LIBRARY_NAME ${PROJECT_NAME}
   SOURCES      "${sources}" )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_memory
   EXPORT draco-targets
   DESTINATION lib )
install( FILES ${headers}  DESTINATION include/memory )

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
