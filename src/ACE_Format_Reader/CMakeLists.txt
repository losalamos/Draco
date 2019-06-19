#-----------------------------*-cmake-*----------------------------------------#
# file   RTT_Format_Reader/CMakeLists.txt
# author Ben Ryan <brryan@lanl.gov>
# date   2019 Jun 19
# brief  Generate build project files for RTT_Format_Reader.
# note   Copyright (C) 2016-2019, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.9.0)
project( ACE_Format_Reader CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB sources *.cc )
file( GLOB headers *.hh )

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library(
   TARGET       Lib_ACE_Format_Reader
   TARGET_DEPS
   LIBRARY_NAME ACE_Format_Reader
   SOURCES      "${sources}"
   HEADERS      "${headers}"
   LIBRARY_NAME_PREFIX "rtt_" )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_ACE_Format_Reader EXPORT draco-targets DESTINATION
  ${DBSCFGDIR}lib )
install( FILES ${headers} DESTINATION ${DBSCFGDIR}include/ACE_Format_Reader )

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
# End ACE_Format_Reader/CMakeLists.txt
#------------------------------------------------------------------------------#
