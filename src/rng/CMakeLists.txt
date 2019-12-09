#-----------------------------*-cmake-*----------------------------------------#
# file   rng/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for rng.
# note   Copyright (C) 2010-2019, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.9.0)
project( rng CXX )

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #
# default to using C++11 features of Random123.
get_target_property( cxxstd Lib_dsxx CXX_STANDARD )
if( cxxstd AND NOT ${cxxstd} STREQUAL 98 )
  set( R123_USE_CXX11 1 )
endif()
configure_file( config.h.in ${PROJECT_BINARY_DIR}/rng/config.h )
set( pkg_config_h "$<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/rng/config.h>")
#\bug can't add "$<INSTALL_INTERFACE:$<CONFIG>/include/rng/config.h>"

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB sources_glob RELATIVE "${PROJECT_SOURCE_DIR}" *.cc *.c )
file( GLOB headers_glob RELATIVE "${PROJECT_SOURCE_DIR}" *.hh *.h *.hpp)

# prefix header files with appropriate paths (build vs install).
foreach( file ${headers_glob} )    
  list(APPEND headers_public "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/${file}>")
  #\bug can't add $<INSTALL_INTERFACE:$<CONFIG>/include/rng>
endforeach()
list(APPEND headers_public "${pkg_config_h}" )
list(APPEND headers_glob ${PROJECT_BINARY_DIR}/rng/config.h ) # install these!

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library(
   TARGET       Lib_rng
   TARGET_DEPS  "Lib_dsxx;GSL::gsl"
   LIBRARY_NAME ${PROJECT_NAME} 
   VENDOR_INCLUDE_DIRS "${RANDOM123_INCLUDE_DIR}" )
target_sources( Lib_rng PUBLIC "${headers_public}" PRIVATE "${sources_glob}" )
set_target_properties(Lib_rng PROPERTIES PUBLIC_HEADER "${headers_public}" )
target_include_directories( Lib_rng PUBLIC 
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}> 
  ${RANDOM123_INCLUDE_DIR} )
  
# \bug I would like to use the following property so that CMake can 
# auto-generate draco-targets.cmake.  However, as of cmake-3.15, I can't write 
# '$<CONFIG>' to exported properties for target Lib_rng.  If this issue fixed 
# by CMake, then all of the machinery in compnent_macros.cmake that sets 
# INTERFACE_INCLUDE_DIRECTORIES in draco-config.cmake can be retired.
# \sa https://gitlab.kitware.com/cmake/cmake/issues/19791
#
#  INTERFACE_INCLUDE_DIRECTORIES $<INSTALL_INTERFACE:$<CONFIG>>/include>"

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_rng EXPORT draco-targets DESTINATION ${DBSCFGDIR}lib )
install( FILES ${headers_glob} DESTINATION ${DBSCFGDIR}include/rng )

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
# End rng/CMakeLists.txt
#------------------------------------------------------------------------------#
