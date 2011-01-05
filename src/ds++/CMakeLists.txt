# $Id$

cmake_minimum_required(VERSION 2.6)
project( dsxx CXX )

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #

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

# Make the header files available in the IDE.
if( MSVC_IDE )
   list( APPEND sources ${headers} )
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

# add_component_library( target_name library_name "$sources" )
add_component_library( Lib_dsxx ds++ "${sources}" )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_dsxx DESTINATION lib )
install( FILES ${headers} ${template_implementations} DESTINATION include/ds++ )

# ---------------------------------------------------------------------------- #
# Unit tests
# ---------------------------------------------------------------------------- #

if( BUILD_TESTING )
   add_subdirectory( test )
endif()   
