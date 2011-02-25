# $Id$

cmake_minimum_required(VERSION 2.6)
project( diagnostics CXX )

##---------------------------------------------------------------------------##
# Special options for this component
##---------------------------------------------------------------------------##

# From autoconf build system: "yes" --> 7, "no" --> 0
set( DRACO_DIAGNOSTICS "0"
   CACHE STRING "Diagnostics verbosity level (0-7). Default is 0." )
# From autoconf build system: "yes" --> 2, "no" --> 0
set( DRACO_TIMING "0"
   CACHE STRING "Timing detail level (0-2). Default is 0." )

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
#list( REMOVE_ITEM headers ${template_implementations} )

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
   ${dsxx_BINARY_DIR}          # ds++/config.h
   ${MPI_INCLUDE_PATH} )

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library( Lib_diagnostics ${PROJECT_NAME} "${sources}" )
add_dependencies( Lib_diagnostics Lib_c4 Lib_dsxx )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_diagnostics DESTINATION lib )
install( FILES ${headers}  DESTINATION include/diagnostics )

# ---------------------------------------------------------------------------- #
# Unit tests
# ---------------------------------------------------------------------------- #

if( BUILD_TESTING )
   add_subdirectory( test )
endif()   
