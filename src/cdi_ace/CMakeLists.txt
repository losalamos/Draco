#-----------------------------*-cmake-*----------------------------------------#
# file   cdi_ace/CMakeLists.txt
# author B. R. Ryan <brryan@lanl.gov>
# date   2019 Jun 19
# brief  Generate build project files for cdi_ace.
# note   Copyright (C) 2016-2019, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.9.0)
project( cdi_ace CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB sources *.cc )
file( GLOB headers *.hh )

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library(
   TARGET       Lib_cdi_ace
   TARGET_DEPS  Lib_dsxx
   LIBRARY_NAME cdi_ace
   SOURCES      "${sources}"
   HEADERS      "${headers}"
   LIBRARY_NAME_PREFIX "rtt_" )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_cdi_ace EXPORT draco-targets DESTINATION
  ${DBSCFGDIR}lib )
install( FILES ${headers} DESTINATION ${DBSCFGDIR}include/cdi_ace )

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
# End cdi_ace/CMakeLists.txt
#------------------------------------------------------------------------------#
