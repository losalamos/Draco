#-----------------------------*-cmake-*----------------------------------------#
# file   cdi_ndi/CMakeLists.txt
# brief  Instructions for building cdi_ndi Makefiles.
# note   Copyright (C) 2020 Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.9.0)
project( cdi_ndi CXX )

if( NDI_FOUND AND NOT DEFINED ENV{TRAVIS} )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB sources *.cc )
file( GLOB headers *.hh )
#list( REMOVE_ITEM sources ${PROJECT_SOURCE_DIR}/QueryNDI.cc )

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library(
   TARGET       Lib_cdi_ndi
   TARGET_DEPS  "Lib_cdi;Lib_rng;NDI::ndi"
   LIBRARY_NAME cdi_ndi
   SOURCES      "${sources}"
   HEADERS      "${headers}" )
target_include_directories( Lib_cdi_ndi PUBLIC ${NDI_INCLUDE_DIR} )

#add_component_executable(
#  TARGET      Exe_QueryNDI
#  TARGET_DEPS Lib_cdi_ndi
#  SOURCES     ${PROJECT_SOURCE_DIR}/QueryNDI.cc
#  PREFIX       Draco )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_cdi_ndi  EXPORT draco-targets DESTINATION
  ${DBSCFGDIR}lib )
#install( TARGETS Exe_QueryNDI EXPORT draco-targets DESTINATION
#  ${DBSCFGDIR}bin )
install( FILES ${headers} DESTINATION ${DBSCFGDIR}include/cdi_ndi )
if( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  install(FILES $<TARGET_PDB_FILE:Lib_cdi_eospac> DESTINATION ${DBSCFGDIR}lib
    OPTIONAL)
endif()

# ---------------------------------------------------------------------------- #
# Unit tests
# ---------------------------------------------------------------------------- #

if( BUILD_TESTING AND NOT DEFINED ENV{APPVEYOR} )
  add_subdirectory( test )
endif()

# ---------------------------------------------------------------------------- #
# Autodoc
# ---------------------------------------------------------------------------- #

process_autodoc_pages()

endif() # NDI_FOUND

# ---------------------------------------------------------------------------- #
# End cdi_ndi/CMakeLists.txt
# ---------------------------------------------------------------------------- #
