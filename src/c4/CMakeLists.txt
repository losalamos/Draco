#-----------------------------*-cmake-*----------------------------------------#
# file   c4/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for c4
# note   Copyright (C) 2016, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.9.0)
project( c4 CXX )

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# many c4 and mpi values are set in config/setupMPI.cmake.
# ---------------------------------------------------------------------------- #
configure_file( config.h.in ${PROJECT_BINARY_DIR}/c4/config.h )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #
file( GLOB template_implementations *.t.hh )
file( GLOB sources *.cc )
file( GLOB explicit_instantiations *_pt.cc )
file( GLOB headers *.hh *.h)
list( APPEND headers ${PROJECT_BINARY_DIR}/c4/config.h )
list( REMOVE_ITEM headers ${template_implementations} )

# Make the header files available in the IDE.
if( MSVC_IDE OR ${CMAKE_GENERATOR} MATCHES Xcode )
   list( APPEND sources ${headers} ${template_implementations} )
endif()

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #
include_directories( ${PROJECT_SOURCE_DIR}       # sources
                     ${PROJECT_BINARY_DIR}       # c4/config.h
                     ${draco_src_dir_SOURCE_DIR} # ds++ header files
                     ${dsxx_BINARY_DIR}          # ds++/config.h
                     )

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #
set( c4_VENDOR_LIBS "${MPI_CXX_LIBRARIES}")
if( PAPI_LIBRARY )
   list( APPEND c4_VENDOR_LIBS ${PAPI_LIBRARY} )
endif()

add_component_library(
   TARGET       Lib_c4
   TARGET_DEPS  Lib_dsxx
   LIBRARY_NAME c4
   SOURCES      "${sources}"
   VENDOR_LIST  "MPI"
   VENDOR_LIBS  "${c4_VENDOR_LIBS}"
   VENDOR_INCLUDE_DIRS "${MPI_CXX_INCLUDE_PATH}"
   )

# xthi binary
if( NOT WIN32 AND NOT APPLE )
  add_subdirectory( bin )
endif()

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #
install( TARGETS Lib_c4 EXPORT draco-targets DESTINATION ${DBSCFGDIR}lib )
install( FILES ${headers} DESTINATION ${DBSCFGDIR}include/c4 )

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
