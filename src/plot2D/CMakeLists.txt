# $Id$

cmake_minimum_required(VERSION 2.6)
project( plot2D CXX )

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #

# Create config.h
if( GRACE_FOUND )
   set( GRACE_H ${Grace_INCLUDE_DIR}/grace_np.h )
endif()
configure_file( config.h.in ${PROJECT_BINARY_DIR}/plot2D/config.h )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

#file( GLOB template_implementations *.t.hh *.i.hh )
file( GLOB sources *.cc )
#file( GLOB explicit_instantiations *_pt.cc )
file( GLOB headers *.hh *.h )
list( APPEND headers ${PROJECT_BINARY_DIR}/plot2D/config.h )
#list( REMOVE_ITEM headers ${template_implementations} )

# Make the header files available in the IDE.
if( MSVC_IDE OR ${CMAKE_GENERATOR} MATCHES Xcode )
   list( APPEND sources ${headers} ${template_implementations} )
endif()

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories( 
   ${PROJECT_SOURCE_DIR}   # sources
   ${PROJECT_BINARY_DIR}   # config.h
   ${draco_src_dir_SOURCE_DIR} # ds++ header files
   ${dsxx_BINARY_DIR}          # ds++/config.h
)

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library( Lib_plot2D plot2D "${sources}" )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_plot2D DESTINATION lib )
install( 
   FILES       ${headers} 
               # ${template_implementations} 
   DESTINATION include/plot2D )

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

