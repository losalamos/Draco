#-----------------------------*-cmake-*----------------------------------------#
# file   config/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for mesh_element/test.
# note   Copyright (C) 2016-2019, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.9.0)
project( mesh_element CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB sources *.cc )
file( GLOB explicit_instantiations *_pt.cc )
file( GLOB headers *.hh )

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library(
   TARGET       Lib_mesh_element
   TARGET_DEPS  Lib_dsxx
   LIBRARY_NAME ${PROJECT_NAME}
   HEADERS      "${headers}"
   SOURCES      "${sources}"  )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_mesh_element EXPORT draco-targets DESTINATION
  ${DBSCFGDIR}lib )
install( FILES ${headers} DESTINATION ${DBSCFGDIR}include/mesh_element )
if( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  install(FILES $<TARGET_PDB_FILE:Lib_mesh_element> DESTINATION ${DBSCFGDIR}lib
    OPTIONAL)
endif()

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
# End mesh_element/CMakeLists.txt
# ---------------------------------------------------------------------------- #
