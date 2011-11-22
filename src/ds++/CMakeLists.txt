# $Id$

cmake_minimum_required(VERSION 2.6)
project( dsxx CXX )

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #

# Store platform information in config.h
if( UNIX )
  set( draco_isLinux 1 )
  set( DRACO_UNAME "Linux" )
elseif( WIN32 )
  set( draco_isWin 1 )
  set( DRACO_UNAME "Windows" )
elseif( OSF1 )
  set( draco_isOSF1 1 )
  set( DRACO_UNAME "OSF1" )
elseif( APPLE )
  set( draco_isDarwin 1)
  set( DRACO_UNAME "Darwin" )
else()
  set( draco_isAIX 1 )
  set( DRACO_UNAME "AIX" )
endif()
# Special setup for catamount
if( ${CMAKE_SYSTEM_NAME} MATCHES "Catamount" )
   set( draco_isLinux_with_aprun 1 )
   set( draco_isCatamount 1 )
   set( DRACO_UNAME "Catamount" )
endif()

# Platform checks for gethostname()
include( CheckFunctionExists)
check_function_exists( gethostname HAVE_GETHOSTNAME )
include( CheckSymbolExists )
check_symbol_exists( HOST_NAME_MAX "unistd.h;limits.h" HAVE_HOST_NAME_MAX )
if( NOT HAVE_HOST_NAME_MAX )
   unset( HAVE_GETHOSTNAME )
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
   ${PROJECT_SOURCE_DIR}/DBC_Ptr.hh )

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

# ---------------------------------------------------------------------------- #
# Autodoc
# ---------------------------------------------------------------------------- #

process_autodoc_pages()

