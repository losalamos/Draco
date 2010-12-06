# $Id$

cmake_minimum_required(VERSION 2.6)
project( special_functions CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB template_implementations *.t.hh *.i.hh )
file( GLOB sources *.cc )
file( GLOB explicit_instantiations *_pt.cc )
file( GLOB headers *.hh )
list( REMOVE_ITEM headers ${template_implementations} )

# Make the header files available in the IDE.
if( MSVC_IDE )
   list( APPEND sources ${headers} )
endif()

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories( 
   ${PROJECT_SOURCE_DIR}       # sources
   ${ode_BINARY_DIR}           # ode/config.h
   ${draco_src_dir_SOURCE_DIR} # ds++ header files
   ${dsxx_BINARY_DIR}          # ds++/config.h
   ${GSL_INCLUDE_DIRS}
   )

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library( Lib_special_functions ${PROJECT_NAME} "${sources}" )
add_dependencies( Lib_special_functions Lib_units Lib_ode Lib_traits)

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_special_functions DESTINATION lib )
install( FILES ${headers} DESTINATION include/special_functions )

# ---------------------------------------------------------------------------- #
# Unit tests
# ---------------------------------------------------------------------------- #

if( DRACO_BUILD_TESTS )
 add_subdirectory( test )
endif()   
