#-----------------------------*-cmake-*----------------------------------------#
# file   src/c4/bin/CMakeLists.txt
# brief  Instructions for building xthi Makefiles
# note   Copyright (C) 2017-2020 Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.5.0)
project( c4_bin CXX )

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_executable(
  TARGET      Exe_xthi
  TARGET_DEPS Lib_c4
  SOURCES     ${PROJECT_SOURCE_DIR}/xthi.cc
  PREFIX       Draco )
#target_include_directories( Exe_xthi
#  PUBLIC $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}> )

add_component_executable(
  TARGET      Exe_ythi
  TARGET_DEPS Lib_c4
  SOURCES     ${PROJECT_SOURCE_DIR}/ythi.cc
  PREFIX       Draco )
#target_include_directories( Exe_ythi
#  PUBLIC $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}> )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Exe_xthi EXPORT draco-targets DESTINATION ${DBSCFGDIR}bin )
install( TARGETS Exe_ythi EXPORT draco-targets DESTINATION ${DBSCFGDIR}bin )

# ---------------------------------------------------------------------------- #
# End cdi_eospac/CMakeLists.txt
# ---------------------------------------------------------------------------- #
