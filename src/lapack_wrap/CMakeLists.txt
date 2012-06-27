# $Id$

cmake_minimum_required(VERSION 2.6)
project( lapack_wrap CXX )

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #
if( NOT lapack_FOUND )
   message( STATUS
      "NOTE: Skipping lapack_wrap component because LAPACK is not found." )
else()

include(FortranCInterface)
FortranCInterface_HEADER( lapack_wrap/config.h MACRO_NAMESPACE "FC_")
# This creates a "lapack_wrap.h" header taht defines mangling macros
# FC_GLOBAL(), FC_GLOBAL_() FC_MODULE(), and FC_MODULE().

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

endif( NOT lapack_FOUND )

# ---------------------------------------------------------------------------- #
# Autodoc
# ---------------------------------------------------------------------------- #

process_autodoc_pages()
