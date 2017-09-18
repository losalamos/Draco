#-----------------------------*-cmake-*----------------------------------------#
# file   src/c4/bin/CMakeLists.txt
# brief  Instructions for building xthi Makefiles
# note   Copyright (C) 2017 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.5.0)
project( c4_bin CXX )

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories(
   ${PROJECT_SOURCE_DIR}      # component headers
   ${PROJECT_BINARY_DIR}      # c4/config.h
   ${PROJECT_SOURCE_DIR}/..   # ds++ and cdi header files
   ${dsxx_BINARY_DIR}         # ds++/config.h
)

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_executable(
  TARGET      Exe_xthi
  TARGET_DEPS Lib_c4
  SOURCES     ${PROJECT_SOURCE_DIR}/xthi.cc
  PREFIX       Draco
  )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Exe_xthi EXPORT draco-targets DESTINATION ${DBSCFGDIR}bin )

# ---------------------------------------------------------------------------- #
# End cdi_eospac/CMakeLists.txt
# ---------------------------------------------------------------------------- #