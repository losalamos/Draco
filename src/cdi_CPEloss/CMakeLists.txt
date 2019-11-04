#-----------------------------*-cmake-*----------------------------------------#
# file   cdi_CPEloss/CMakeLists.txt
# brief  Instructions for building cdi_eospac Makefiles.
# note   Copyright (C) 2016-2019 Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.9.0)
project( cdi_CPEloss CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB sources *.cc )
file( GLOB headers *.hh )

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library(
   TARGET       Lib_cdi_CPEloss
   TARGET_DEPS  "Lib_cdi"
   LIBRARY_NAME ${PROJECT_NAME}
   SOURCES      "${sources}"
   HEADERS      "${headers}" )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_cdi_CPEloss  EXPORT draco-targets DESTINATION
  ${DBSCFGDIR}lib )
install( FILES ${headers} DESTINATION ${DBSCFGDIR}include/cdi_CPEloss )

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

# ---------------------------------------------------------------------------- #
# End cdi_CPEloss/CMakeLists.txt
# ---------------------------------------------------------------------------- #
