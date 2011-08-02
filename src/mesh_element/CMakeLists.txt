# $Id$

cmake_minimum_required(VERSION 2.6)
project( mesh_element CXX )

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #

# configure_file( config.h.in ${PROJECT_BINARY_DIR}/mesh_element/config.h )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

# file( GLOB template_implementations *.t.hh *.i.hh )
file( GLOB sources *.cc )
file( GLOB explicit_instantiations *_pt.cc )
file( GLOB headers *.hh )
#list( REMOVE_ITEM headers ${template_implementations} )

# Make the header files available in the IDE.
if( MSVC_IDE OR ${CMAKE_GENERATOR} MATCHES Xcode )
   list( APPEND sources ${headers} )
endif()

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories( 
   ${PROJECT_SOURCE_DIR}      # sources
   ${PROJECT_SOURCE_DIR}/..   # ds++ header files
   ${dsxx_BINARY_DIR}         # ds++/config.h
) 

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library( Lib_mesh_element ${PROJECT_NAME} "${sources}"  )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_mesh_element DESTINATION lib )
install( FILES ${headers} DESTINATION include/mesh_element )

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


