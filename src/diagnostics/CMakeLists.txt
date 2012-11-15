#-----------------------------*-cmake-*----------------------------------------#
# file   config/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for diagnostics.
# note   Copyright (C) 2010-2012, Los Alamos National Security
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 2.6)
project( diagnostics CXX )

##---------------------------------------------------------------------------##
# Special options for this component
##---------------------------------------------------------------------------##

# From autoconf build system: "yes" --> 7, "no" --> 0
set( DRACO_DIAGNOSTICS "0"
   CACHE STRING "Diagnostics verbosity level (0-7). Default is 0." )
# provide a constrained drop down menu for cmake-gui
set_property( CACHE DRACO_DIAGNOSTICS PROPERTY STRINGS 0 1 2 3 4 5 6 7 )
# From autoconf build system: "yes" --> 2, "no" --> 0
set( DRACO_TIMING "0"
   CACHE STRING "Timing detail level (0-2). Default is 0." )
# provide a constrained drop down menu for cmake-gui
set_property( CACHE DRACO_TIMING PROPERTY STRINGS 0 1 2 )

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #

configure_file( config.h.in ${PROJECT_BINARY_DIR}/diagnostics/config.h )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB template_implementations *.t.hh *.i.hh )
file( GLOB sources *.cc )
file( GLOB explicit_instantiations *_pt.cc )
file( GLOB headers *.hh )
list( APPEND headers ${PROJECT_BINARY_DIR}/diagnostics/config.h )
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
   ${c4_BINARY_DIR}   # c4/config.h
   )

# setup binary directory
file( MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/bin )

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library( 
   TARGET       Lib_diagnostics
   TARGET_DEPS  Lib_c4
   LIBRARY_NAME ${PROJECT_NAME}
   SOURCES      "${sources}" )

add_executable( Exe_draco_info ${PROJECT_SOURCE_DIR}/draco_info.cc )
target_link_libraries( Exe_draco_info Lib_dsxx  )
set_target_properties( Exe_draco_info
   PROPERTIES OUTPUT_NAME bin/draco_info )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_diagnostics  
   EXPORT draco-targets
   DESTINATION lib )
install( TARGETS Exe_draco_info DESTINATION bin )
install( FILES ${headers}  DESTINATION include/diagnostics )

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
