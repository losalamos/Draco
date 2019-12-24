#-----------------------------*-cmake-*----------------------------------------#
# file   cdi_eospac/CMakeLists.txt
# brief  Instructions for building cdi_eospac Makefiles.
# note   Copyright (C) 2016-2019 Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.9.0)
project( cdi_eospac CXX )

if( EOSPAC_FOUND )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB sources *.cc )
file( GLOB headers *.hh )
list( REMOVE_ITEM sources ${PROJECT_SOURCE_DIR}/QueryEospac.cc )

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library(
   TARGET       Lib_cdi_eospac
   TARGET_DEPS  "Lib_cdi;EOSPAC::eospac"
   LIBRARY_NAME cdi_eospac
   SOURCES      "${sources}"
   HEADERS      "${headers}" )
target_include_directories( Lib_cdi_eospac PUBLIC ${EOSPAC_INCLUDE_DIR} )

add_component_executable(
  TARGET      Exe_QueryEospac
  TARGET_DEPS Lib_cdi_eospac
  SOURCES     ${PROJECT_SOURCE_DIR}/QueryEospac.cc
  PREFIX       Draco )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_cdi_eospac  EXPORT draco-targets DESTINATION
  ${DBSCFGDIR}lib )
install( TARGETS Exe_QueryEospac EXPORT draco-targets DESTINATION
  ${DBSCFGDIR}bin )
install( FILES ${headers} DESTINATION ${DBSCFGDIR}include/cdi_eospac )
if( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  install(FILES $<TARGET_PDB_FILE:Lib_cdi_eospac> DESTINATION ${DBSCFGDIR}lib
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

endif() # EOSPAC_FOUND

# ---------------------------------------------------------------------------- #
# End cdi_eospac/CMakeLists.txt
# ---------------------------------------------------------------------------- #
