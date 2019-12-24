#-----------------------------*-cmake-*----------------------------------------#
# file   diagnostics/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for diagnostics.
# note   Copyright (C) 2016-2019, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.9.0)
project( diagnostics CXX )

##---------------------------------------------------------------------------##
# Special options for this component
##---------------------------------------------------------------------------##

# DRACO_DIAGNOSICS is set in ds++ (ds++/config.h).

# From autoconf build system: "yes" --> 2, "no" --> 0
set( DRACO_TIMING "0"
   CACHE STRING "Timing detail level (0-2). Default is 0." )
# provide a constrained drop down menu for cmake-gui
set_property( CACHE DRACO_TIMING PROPERTY STRINGS 0 1 2 )

# Are the memory features of getrusage to be activated?
option( USE_PROCMON "Use PROCMON process (cpu/memory) monitoring." OFF )
# If USE_PROCMON is requested but /proc/meminfo is unavailable, abort
if( USE_PROCMON )
   if( NOT EXISTS "/proc/meminfo" )
      message( FATAL_ERROR "USE_PROCMON=ON, but /proc/meminfo unavailable." )
   endif()
endif()
if( NOT ${HAVE_SYS_TYPES} )
   message( FATAL_ERROR "diagnostics/Procmon.hh requires sys/types.h.")
endif()

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #

configure_file( config.h.in ${PROJECT_BINARY_DIR}/diagnostics/config.h )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB sources *.cc )
file( GLOB headers *.hh )
list( APPEND headers ${PROJECT_BINARY_DIR}/diagnostics/config.h )

# Make the header files available in the IDE.
if( CMAKE_CONFIGURATION_TYPES )
  add_definitions(-DCMAKE_BUILD_TYPE=\"${CMAKE_CFG_INTDIR}\")
endif()

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library(
  TARGET       Lib_diagnostics
  TARGET_DEPS  "Lib_c4"
  LIBRARY_NAME ${PROJECT_NAME}
  SOURCES      "${sources}"
  HEADERS      "${headers}"
  PREFIX       Draco
  )
target_include_directories( Lib_diagnostics
  PUBLIC $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}> )

add_component_executable(
  TARGET      Exe_draco_info
  TARGET_DEPS Lib_diagnostics
  EXE_NAME    draco_info
  SOURCES     ${PROJECT_SOURCE_DIR}/draco_info_main.cc
  HEADERS     ${PROJECT_SOURCE_DIR}/draco_info.hh
  PREFIX       Draco
  )

if( QT_FOUND )
  add_subdirectory( qt )
endif()

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_diagnostics EXPORT draco-targets DESTINATION
  ${DBSCFGDIR}lib )
install( TARGETS Exe_draco_info  EXPORT draco-targets DESTINATION
  ${DBSCFGDIR}bin )
install( FILES ${headers}  DESTINATION ${DBSCFGDIR}include/diagnostics )
if( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  install(FILES $<TARGET_PDB_FILE:Lib_diagnostics> DESTINATION ${DBSCFGDIR}lib
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

#------------------------------------------------------------------------------#
# End diagnostics/CMakeLists.txt
#------------------------------------------------------------------------------#
