# $Id$

cmake_minimum_required(VERSION 2.6)
project( cdi_gandolf CXX )

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #

if( NOT GANDOLF_FOUND )
   set( rtt_cdi_gandolf_stub 1 )
endif()

configure_file( config.h.in ${PROJECT_BINARY_DIR}/cdi_gandolf/config.h )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

# file( GLOB template_implementations *.t.hh *.i.hh )
file( GLOB sources *.cc )
file( GLOB explicit_instantiations *_pt.cc )
file( GLOB headers *.hh )
list( APPEND headers ${PROJECT_BINARY_DIR}/cdi_gandolf/config.h )

# Make the header files available in the IDE.
if( MSVC_IDE )
   list( APPEND sources ${headers} )
endif()

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories( 
   ${PROJECT_SOURCE_DIR}      # sources
   ${PROJECT_SOURCE_DIR}/..   # header files
   ${PROJECT_BINARY_DIR}      # cdi_gandolf/config.h
   ${dsxx_BINARY_DIR}         # ds++/config.h
) 

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_library( Lib_cdi_gandolf ${DRACO_LIBRARY_TYPE} ${sources}  )
add_dependencies( Lib_cdi_gandolf Lib_cdi Lib_dsxx )
if( "${DRACO_LIBRARY_TYPE}" MATCHES "SHARED" )
   set_target_properties( Lib_cdi_gandolf 
      PROPERTIES 
      # Provide compile define macro to enable declspec(dllexport) linkage.
      COMPILE_DEFINITIONS BUILDING_DLL 
      # Use custom library naming
      OUTPUT_NAME rtt_cdi_gandolf
      )
else()
   set_target_properties( Lib_cdi_gandolf 
      PROPERTIES 
      # Use custom library naming
      OUTPUT_NAME rtt_cdi_gandolf
      )
endif()

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_cdi_gandolf DESTINATION lib )
install( FILES ${headers} DESTINATION include/cdi_gandolf )

# ---------------------------------------------------------------------------- #
# Unit tests
# ---------------------------------------------------------------------------- #

if( DRACO_BUILD_TESTS AND EXISTS ${PROJECT_SOURCE_DIR}/test/CMakeLists.txt )
   add_subdirectory( test )
endif()   

