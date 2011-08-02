# $Id$

cmake_minimum_required(VERSION 2.6)
project( RTT_Format_Reader CXX )

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

include_directories( ${PROJECT_SOURCE_DIR}       # sources
                     ${draco_src_dir_SOURCE_DIR} # ds++ header files
                     ${dsxx_BINARY_DIR}          # ds++/config.h
                     ${MPI_INCLUDE_PATH}
)

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library( Lib_RTT_Format_Reader RTT_Format_Reader "${sources}" )
add_dependencies( Lib_RTT_Format_Reader
   Lib_meshReaders
   Lib_mesh_element 
   Lib_dsxx )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_RTT_Format_Reader DESTINATION lib )
install( FILES ${headers} DESTINATION include/RTT_Format_Reader )

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


