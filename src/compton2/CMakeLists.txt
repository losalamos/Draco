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

set( sources
    #${PROJECT_SOURCE_DIR}/Sparse_Compton_Matrix.cc
    ${PROJECT_SOURCE_DIR}/Compton2.cc
)

set( headers
    #${PROJECT_SOURCE_DIR}/Sparse_Compton_Matrix.hh
  ${PROJECT_SOURCE_DIR}/Compton2.hh
)

# -------------------------------------------------------------------------- #
# Build package library
# -------------------------------------------------------------------------- #

add_component_library(
  TARGET       Lib_compton2
  TARGET_DEPS  "Lib_c4;Lib_units;Lib_cdi_analytic"
  INCLUDE_DIRS "PUBLIC;$<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>"
  LIBRARY_NAME compton2
  SOURCES      "${sources}"
  HEADERS      "${headers}" )

add_component_executable(
  TARGET      Exe_CSKrw
  TARGET_DEPS "Lib_units;Lib_compton2"
  SOURCES     ${PROJECT_SOURCE_DIR}/cskrw.cc)

# Copy necessary dll files to the build directory
if( MSVC AND DRACO_LIBRARY_TYPE STREQUAL SHARED )
add_custom_command( TARGET Lib_compton2 POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
  $<TARGET_FILE_DIR:Lib_compton2> )
endif()

# -------------------------------------------------------------------------- #
# Installation instructions
# -------------------------------------------------------------------------- #
set( CMAKE_INSTALL_DO_STRIP NO ) # ??
install( FILES ${headers} DESTINATION ${DBSCFGDIR}include/compton2 )
install( TARGETS Exe_CSKrw EXPORT draco-targets DESTINATION ${DBSCFGDIR}bin )
install( FILES python/common_compton.py python/csk_reader.py python/merge_csk.py DESTINATION ${DBSCFGDIR}bin )

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
