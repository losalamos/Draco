#-----------------------------*-cmake-*----------------------------------------#
# file   rng/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for rng.
# note   Copyright (C) 2010-2019, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.9.0)
project( rng CXX )

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #
# default to using C++11 features of Random123.
set( R123_USE_CXX11 1 )
configure_file( config.h.in ${PROJECT_BINARY_DIR}/rng/config.h )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB sources *.cc *.c )
file( GLOB headers *.hh *.h *.hpp)
list( APPEND headers ${PROJECT_BINARY_DIR}/rng/config.h )

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library(
   TARGET       Lib_rng
   TARGET_DEPS  "Lib_dsxx;GSL::gsl"
   LIBRARY_NAME ${PROJECT_NAME}
   SOURCES      "${sources}"
   HEADERS      "${headers}"
   )
target_include_directories( Lib_rng
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}> # config.h
    $<BUILD_INTERFACE:${RANDOM123_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:${RANDOM123_INCLUDE_DIR}> )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_rng EXPORT draco-targets DESTINATION ${DBSCFGDIR}lib )
install( FILES ${headers} DESTINATION ${DBSCFGDIR}include/rng )

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
# End rng/CMakeLists.txt
#------------------------------------------------------------------------------#
