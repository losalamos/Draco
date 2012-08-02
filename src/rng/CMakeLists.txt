# $Id$

cmake_minimum_required(VERSION 2.6)
project( rng CXX )

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #

option( ENABLE_RNG_NR "Enable non reproducable random number generator" OFF )
if( ENABLE_RNG_NR )
   set( RNG_NR 1 )
endif()
configure_file( config.h.in ${PROJECT_BINARY_DIR}/rng/config.h )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB sources *.cc *.c )
file( GLOB headers *.hh *.h )
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
   ${GSL_INCLUDE_DIR}
) 

if( RANDOM123_FOUND )
   include_directories( ${RANDOM123_INCLUDE_DIR} )
endif()

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library( Lib_rng ${PROJECT_NAME} "${sources}"  )
target_link_libraries( Lib_rng ${GSL_LIBRARIES}  )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_rng DESTINATION lib )
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


