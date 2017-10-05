#-----------------------------*-cmake-*----------------------------------------#
# file   src/compton/CMakeLists.txt
# author Kendra Keady
# date   2017 February 28
# brief  Instructions for building compton Makefiles.
# note   Copyright (C) 2017 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.9.0)
project( compton CXX )

# This packages requires libcompton.  See draco/config/vendor_libraries.cmake
# and draco/config/FindCOMPTON.cmake to manage how the build system discovers
# this TPL.  If not found, only a stub header file will be installed.

configure_file( ${PROJECT_SOURCE_DIR}/config.h.in
  ${PROJECT_BINARY_DIR}/config.h )

set( headers
  ${PROJECT_BINARY_DIR}/config.h
  ${PROJECT_SOURCE_DIR}/Compton.hh )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

if( COMPTON_FOUND )

set( sources ${PROJECT_SOURCE_DIR}/Compton.cc )

  # Make the header files available in the IDE.
  if( MSVC_IDE OR ${CMAKE_GENERATOR} MATCHES Xcode )
    list( APPEND sources ${headers} )
  endif()

  # -------------------------------------------------------------------------- #
  # Directories to search for include directives
  # -------------------------------------------------------------------------- #

  include_directories(
    ${Draco_SOURCE_DIR}/src    # compton/Compton.hh ds++/Assert.hh
    ${Draco_BINARY_DIR}/src    # compton/config.h
    ${dsxx_BINARY_DIR}         # ds++/config.h
    ${c4_BINARY_DIR}           # c4/config.h
    )

  # -------------------------------------------------------------------------- #
  # Build package library
  # -------------------------------------------------------------------------- #

  add_component_library(
    TARGET       Lib_compton
    TARGET_DEPS  "Lib_c4;Lib_dsxx;COMPTON::compton"
    LIBRARY_NAME compton
    SOURCES      "${sources}"
)

  install( TARGETS Lib_compton  EXPORT draco-targets DESTINATION
    ${DBSCFGDIR}lib )

endif() # COMPTON_FOUND

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( FILES ${headers} DESTINATION ${DBSCFGDIR}include/compton )

# ---------------------------------------------------------------------------- #
# Unit tests
# ---------------------------------------------------------------------------- #

if( BUILD_TESTING AND COMPTON_FOUND )
  add_subdirectory( test )
endif()

# ---------------------------------------------------------------------------- #
# Autodoc
# ---------------------------------------------------------------------------- #

process_autodoc_pages()

#------------------------------------------------------------------------------#
# End compton/CMakeLists.txt
#------------------------------------------------------------------------------#
