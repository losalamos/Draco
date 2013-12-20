#-----------------------------*-cmake-*----------------------------------------#
# file   config/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for rng.
# note   Copyright (C) 2010-2012, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#
cmake_minimum_required( VERSION 2.8 )
project( rng CXX )

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #

option( ENABLE_RNG_NR "Enable non reproducable random number generator" OFF )
if( ENABLE_RNG_NR )
   set( RNG_NR 1 )
endif()
if( HAS_CXX11_EXPLICIT_CONVERSION AND HAS_CXX11_STATIC_ASSERT AND
      HAS_CXX11_TYPE_TRAITS AND HAS_CXX11_UNRESTRICTED_UNIONS AND
      HAS_CXX11_LONG_LONG AND NOT USE_CUDA )
   set( R123_USE_CXX11 1 )
endif()
configure_file( config.h.in ${PROJECT_BINARY_DIR}/rng/config.h )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB sources *.cc *.c )
file( GLOB headers *.hh *.h *.hpp)
list( APPEND headers ${PROJECT_BINARY_DIR}/rng/config.h )

# Make the header files available in the IDE.
if( MSVC_IDE OR ${CMAKE_GENERATOR} MATCHES Xcode )
   list( APPEND sources ${headers} )
endif()

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories( 
   ${PROJECT_SOURCE_DIR}      # sources
   ${PROJECT_BINARY_DIR}      # rng/config.h
   ${PROJECT_SOURCE_DIR}/..   # ds++ header files
   ${dsxx_BINARY_DIR}         # ds++/config.h
) 

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

if( RANDOM123_FOUND )
   set( rng_vendor_include_dirs
      "${GSL_INCLUDE_DIRS};${RANDOM123_INCLUDE_DIR}" )
else()
   set( rng_vendor_include_dirs "${GSL_INCLUDE_DIRS}" )
endif()

add_component_library( 
   TARGET       Lib_rng
   TARGET_DEPS  Lib_dsxx
   LIBRARY_NAME ${PROJECT_NAME}
   SOURCES      "${sources}"  
   VENDOR_LIST         "GSL"
   VENDOR_LIBS         "${GSL_LIBRARIES}"
   VENDOR_INCLUDE_DIRS "${rng_vendor_include_dirs}"
)

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_rng EXPORT draco-targets DESTINATION lib )
install( FILES ${headers} DESTINATION include/rng )

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
