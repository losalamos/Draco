# $Id$

cmake_minimum_required(VERSION 2.6)
project( fit CXX )

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #

configure_file( config.h.in ${PROJECT_BINARY_DIR}/fit/config.h )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB template_implementations *.t.hh *.i.hh )
file( GLOB sources *.cc )
file( GLOB headers *.hh )
list( APPEND headers ${PROJECT_BINARY_DIR}/fit/config.h )
list( REMOVE_ITEM headers ${template_implementations} )

# Make the header files available in the IDE.
if( MSVC_IDE )
   list( APPEND sources ${headers} )
endif()

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories( 
   ${PROJECT_SOURCE_DIR}      # sources
   ${PROJECT_BINARY_DIR}      # fit/config.h
   ${PROJECT_SOURCE_DIR}/..   # ds++ header files
   ${dsxx_BINARY_DIR}         # ds++/config.h
) 

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_library( Lib_fit ${DRACO_LIBRARY_TYPE} ${sources} )
add_dependencies( Lib_fit Lib_linear )
if( "${DRACO_LIBRARY_TYPE}" MATCHES "SHARED" )
   set_target_properties( Lib_fit 
      PROPERTIES 
      # Provide compile define macro to enable declspec(dllexport) linkage.
      COMPILE_DEFINITIONS BUILDING_DLL 
      # Use custom library naming
      OUTPUT_NAME rtt_fit
      )
else()
   set_target_properties( Lib_fit 
      PROPERTIES 
      # Use custom library naming
      OUTPUT_NAME rtt_fit
      )
endif()
#target_link_libraries( Lib_fit Lib_c4 Lib_dsxx )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_fit DESTINATION lib )
install( FILES ${headers} ${template_implementations} DESTINATION include/fit )

# ---------------------------------------------------------------------------- #
# Unit tests
# ---------------------------------------------------------------------------- #

if( DRACO_BUILD_TESTS )
   add_subdirectory( test )
endif()   
