#-----------------------------*-cmake-*----------------------------------------#
# file   src/compton/CMakeLists.txt
# author Andrew Till
# date   2020 May 11
# brief  Instructions for building Compton Makefiles.
# note   Copyright (C) 2017-2020 Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.9.0)
project( compton2 CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( sources "${PROJECT_SOURCE_DIR}/Compton2.cc" )
set( headers
  ${PROJECT_SOURCE_DIR}/Compton2.hh )

# -------------------------------------------------------------------------- #
# Build package library
# -------------------------------------------------------------------------- #

add_component_library(
TARGET       Lib_compton2
TARGET_DEPS  "Lib_c4"
LIBRARY_NAME compton2
SOURCES      "${sources}"
HEADERS      "${headers}" )
# generated include directive files (config.h)
target_include_directories( Lib_compton2
  PUBLIC $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}> )

# Copy necessary dll files to the build directory
if( MSVC AND DRACO_LIBRARY_TYPE STREQUAL SHARED )
add_custom_command( TARGET Lib_compton2 POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different 
  $<TARGET_FILE_DIR:Lib_compton2> )
endif()

# -------------------------------------------------------------------------- #
# Installation instructions
# -------------------------------------------------------------------------- #

install( TARGETS Lib_compton2 EXPORT draco-targets DESTINATION
  ${DBSCFGDIR}lib )
install( FILES ${headers} DESTINATION ${DBSCFGDIR}include/compton2 )
if( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
 install(FILES $<TARGET_PDB_FILE:Lib_compton2> DESTINATION ${DBSCFGDIR}lib
   OPTIONAL)
endif()

# -------------------------------------------------------------------------- #
# Unit tests
# -------------------------------------------------------------------------- #

if( BUILD_TESTING )
  add_subdirectory( test )
endif()

# -------------------------------------------------------------------------- #
# Autodoc
# -------------------------------------------------------------------------- #

process_autodoc_pages()

#------------------------------------------------------------------------------#
# End compton2/CMakeLists.txt
#------------------------------------------------------------------------------#
