#-----------------------------*-cmake-*----------------------------------------#
# file   cdi_ndi/CMakeLists.txt
# brief  Instructions for building cdi_ndi Makefiles.
# note   Copyright (C) 2020 Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.9.0)
project( cdi_ndi CXX )

if( NDI_FOUND AND NOT DEFINED ENV{TRAVIS} )

if( NOT DEFINED NDI_ROOT_DIR AND DEFINED $ENV{NDI_ROOT_DIR})
  set(NDI_ROOT_DIR $ENV{NDI_ROOT_DIR})
endif()
if( NOT DEFINED NDI_ROOT_DIR )
  message(FATAL_ERROR "NDI_ROOT_DIR not set. cdi_ndi requires the variable \
          NDI_ROOT_DIR to be set in the environment or with the -DNDI_ROOT_DIR \
          flag if the ndi library is found.")
endif()

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #
configure_file( config.h.in ${PROJECT_BINARY_DIR}/cdi_ndi/config.h )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB sources *.cc )
file( GLOB headers *.hh )
list( APPEND headers ${PROJECT_BINARY_DIR}/cdi_ndi/config.h )

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library(
   TARGET       Lib_cdi_ndi
   TARGET_DEPS  "Lib_cdi;Lib_rng;Lib_units;NDI::ndi"
   LIBRARY_NAME cdi_ndi
   SOURCES      "${sources}"
   HEADERS      "${headers}" )
target_include_directories( Lib_cdi_ndi PUBLIC ${NDI_INCLUDE_DIR}
   PUBLIC $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}> )

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
