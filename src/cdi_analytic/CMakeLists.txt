#-----------------------------*-cmake-*----------------------------------------#
# file   config/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for cdi_analytic.
# note   Copyright (C) 2016, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.9.0)
project( cdi_analytic CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB sources *.cc )
file( GLOB headers *.hh )

# Make the header files available in the IDE.
if( MSVC_IDE OR ${CMAKE_GENERATOR} MATCHES Xcode )
   list( APPEND sources ${headers} )
endif()

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories(
   ${PROJECT_SOURCE_DIR}      # component headers
   ${PROJECT_SOURCE_DIR}/..   # ds++ and cdi header files
   ${dsxx_BINARY_DIR}         # ds++/config.h
   ${c4_BINARY_DIR}           # c4/config.h
   ${ode_BINARY_DIR}          # ode/config.h
   ${parser_BINARY_DIR}
   ${MPI_CXX_INCLUDE_PATH}
)

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library(
   TARGET       Lib_cdi_analytic
   TARGET_DEPS  "Lib_parser;Lib_roots;Lib_cdi"
   #TARGET_DEPS  "Lib_parser;Lib_roots;Lib_cdi;Lib_ode"
   LIBRARY_NAME cdi_analytic
   SOURCES      "${sources}"  )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_cdi_analytic EXPORT draco-targets DESTINATION ${DBSCFGDIR}lib )
install( FILES ${headers} DESTINATION ${DBSCFGDIR}include/cdi_analytic )

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
