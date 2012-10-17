#-----------------------------*-cmake-*----------------------------------------#
# file   config/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for lapack_wrap.
# note   Copyright (C) 2010-2012, Los Alamos National Security
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 2.6)
project( lapack_wrap CXX )

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #
if( NOT lapack_FOUND )
   message( STATUS
      "NOTE: Skipping lapack_wrap component because LAPACK is not found." )
else()

# cmake --help-module FortranCInterface
# This module automatically detects the API by which C and Fortran
# languages interact.
include(FortranCInterface)

# A function is provided to generate a C header file containing macros
# to mangle symbol names: 
#
# FortranCInterface_HEADER(<file>
#                          [MACRO_NAMESPACE <macro-ns>]
#                          [SYMBOL_NAMESPACE <ns>]
#                          [SYMBOLS [<module>:]<function> ...])
FortranCInterface_HEADER( lapack_wrap/config.h MACRO_NAMESPACE "FC_")
# This creates a "lapack_wrap.h" header that defines mangling macros
# FC_GLOBAL(), FC_GLOBAL_() FC_MODULE(), and FC_MODULE().

# Omit this package is the mangling was not found by
# FortranCInterface.cmake 
if( FortranCInterface_GLOBAL_FOUND AND FortranCInterface_MODULE_FOUND )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB sources *.cc )
file( GLOB headers *.hh )
list( APPEND headers ${PROJECT_BINARY_DIR}/lapack_wrap/config.h )

# Make the header files available in the IDE.
if( MSVC_IDE OR ${CMAKE_GENERATOR} MATCHES Xcode )
   list( APPEND sources ${headers} )
endif()

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories( 
   ${PROJECT_SOURCE_DIR}     # sources
   ${PROJECT_BINARY_DIR}     # lapack_wrap/config.h
   ${PROJECT_SOURCE_DIR}/..  # ds++ header files
   ${dsxx_BINARY_DIR}        # ds++/config.h
) 

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

# No sources!
# add_component_library( Lib_lapack_wrap ${PROJECT_NAME} "${sources}"  )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( FILES ${headers} DESTINATION include/lapack_wrap )

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

# --------------------------------------------------------------------------- #
# FortranCInterface_GLOBAL_FOUND AND FortranCInterface_MODULE_FOUND 
else()
   message("   lapack_wrap - skipped because FortranCInterface failed.")
endif()

# --------------------------------------------------------------------------- #
# lapack_FOUND?
else()
   message("   lapack_wrap - skipped because lapack was not found.")
endif( NOT lapack_FOUND )