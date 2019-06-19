#-----------------------------*-cmake-*----------------------------------------#
# file   nuclear_data/CMakeLists.txt
# author B. R. Ryan <brryan@lanl.gov>
# date   2019 Jun 19
# brief  Generate build project files for nuclear_data.
# note   Copyright (C) 2016-2019, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.9.0)
project( nuclear_data CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB sources *.cc )
file( GLOB headers *.hh )

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library(
   TARGET       Lib_nuclear_data
   TARGET_DEPS
   LIBRARY_NAME nuclear_data
   SOURCES      "${sources}"
   HEADERS      "${headers}"
   LIBRARY_NAME_PREFIX "rtt_" )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_nuclear_data EXPORT draco-targets DESTINATION
  ${DBSCFGDIR}lib )
install( FILES ${headers} DESTINATION ${DBSCFGDIR}include/nuclear_data )

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
# End nuclear_data/CMakeLists.txt
#------------------------------------------------------------------------------#
