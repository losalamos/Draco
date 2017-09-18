#-----------------------------*-cmake-*----------------------------------------#
# file   config/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for rng.
# note   Copyright (C) 2010-2016, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.9.0)
project( rng CXX )

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #
# default to using C++11 features of Random123.
set( R123_USE_CXX11 1 )
# if the current compiler doesn't support some c++11 features, disable
# Random123's use of these extensions.
set( required_features "cxx_explicit_conversions;cxx_static_assert;cxx_unrestricted_unions;cxx_long_long_type;cxx_constexpr;cxx_type_traits")
foreach( feature ${required_features} )
   list( FIND CXX11_FEATURE_LIST ${feature} featurefound )
   if( ${featurefound} STREQUAL "-1" )
      unset(R123_USE_CXX11)
   endif()
endforeach()

# When also linking to CUDA, C++11 doesn't work.
# MSVC 2013 RTM doesn't provide constexpr, even though cmake thinks it does.
if( USE_CUDA OR WIN32 )
   unset(R123_USE_CXX11)
endif()

# At least one report that Intel/15 doesn't work correctly with Random123 with
# c++11 features active.
# if( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel" AND
    # ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 16.0 )
   # unset(R123_USE_CXX11)
# endif()

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

add_component_library(
   TARGET       Lib_rng
   TARGET_DEPS  "Lib_dsxx;GSL::gsl"
   LIBRARY_NAME ${PROJECT_NAME}
   SOURCES      "${sources}"
   VENDOR_LIST  "Random123"
   VENDOR_INCLUDE_DIRS "${RANDOM123_INCLUDE_DIR}"
)

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_rng EXPORT draco-targets DESTINATION ${DBSCFGDIR}lib )
install( FILES ${headers} DESTINATION ${DBSCFGDIR}include/rng )

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
