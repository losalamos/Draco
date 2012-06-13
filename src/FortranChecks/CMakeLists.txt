# $Id: CMakeLists.txt 6475 2012-03-12 22:54:19Z kellyt $

cmake_minimum_required(VERSION 2.6)
project( FortranChecks CXX Fortran )

#include(FortranCInterface)
#FortranCInterface_HEADER( FortranChecks/config.h MACRO_NAMESPACE "FC_")
# This creates a "config.h" header that defines mangling macros
# FC_GLOBAL(), FC_GLOBAL_() FC_MODULE(), and FC_MODULE().

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #

# Review config/platform_checks.cmake for many of these settings.

# Create ds++/config.h
# configure_file( config.h.in ${PROJECT_BINARY_DIR}/ds++/config.h )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB f90sources *.f90 )

#file( GLOB template_implementations *.t.hh *.i.hh )
#file( GLOB sources *.cc )
#file( GLOB explicit_instantiations *_pt.cc )
#file( GLOB headers *.hh )
#list( APPEND headers ${PROJECT_BINARY_DIR}/ds++/config.h )
#list( REMOVE_ITEM headers ${template_implementations} )

# Make the header files available in the IDE.
#if( MSVC_IDE OR ${CMAKE_GENERATOR} MATCHES Xcode )
#   list( APPEND sources ${headers} ${template_implementations} )
#endif()

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

#include_directories( 
#   ${PROJECT_SOURCE_DIR}   # sources
#)
#   ${PROJECT_BINARY_DIR} ) # config.h

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

# add_component_library( target_name library_name "$sources" )
add_component_library( Lib_FortranChecks FortranCheck "${f90sources}" )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

# install( TARGETS Lib_dsxx DESTINATION lib )
#target_link_libraries( Lib_dsxx ${Lib_win_winsock} )
#install( FILES ${headers} ${template_implementations} DESTINATION include/ds++ )

# ---------------------------------------------------------------------------- #
# Unit tests
# ---------------------------------------------------------------------------- #

if( BUILD_TESTING )
   add_subdirectory( test )
   add_subdirectory( ftest )
endif()   

# ---------------------------------------------------------------------------- #
# Autodoc
# ---------------------------------------------------------------------------- #

process_autodoc_pages()

