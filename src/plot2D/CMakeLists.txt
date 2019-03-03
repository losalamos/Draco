#-----------------------------*-cmake-*----------------------------------------#
# file   plot2D/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for plot2D.
# note   Copyright (C) 2016-2019, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.9.0)
project( plot2D CXX )

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #

# Create config.h
if( Grace_FOUND )
   set( GRACE_H ${Grace_INCLUDE_DIR}/grace_np.h )
endif()
configure_file( config.h.in ${PROJECT_BINARY_DIR}/plot2D/config.h )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB sources *.cc )
file( GLOB headers *.hh *.h )
list( APPEND headers ${PROJECT_BINARY_DIR}/plot2D/config.h )

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library(
   TARGET       Lib_plot2D
   TARGET_DEPS  Lib_dsxx
   LIBRARY_NAME plot2D
   SOURCES      "${sources}"
   HEADERS      "${headers}" )
 target_include_directories( Lib_plot2D
   PUBLIC
     ${Grace_INCLUDE_DIR}
     $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}> )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_plot2D EXPORT draco-targets DESTINATION ${DBSCFGDIR}lib )
install( FILES ${headers} DESTINATION ${DBSCFGDIR}include/plot2D )

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

#------------------------------------------------------------------------------#
# End plot2D/CMakeLists.txt
#------------------------------------------------------------------------------#
