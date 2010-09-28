# $Id$

cmake_minimum_required(VERSION 2.6)
project( dsxx CXX )

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #

configure_file( config.h.in ${PROJECT_BINARY_DIR}/ds++/config.h )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB template_implementations *.t.hh *.i.hh )
file( GLOB sources *.cc )
file( GLOB explicit_instantiations *_pt.cc )
file( GLOB headers *.hh )
list( APPEND headers ${PROJECT_BINARY_DIR}/ds++/config.h )
list( REMOVE_ITEM headers ${template_implementations} )

# Make the header files available in the IDE.
if( MSVC_IDE )
   list( APPEND sources ${headers} )
endif()

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories( ${PROJECT_SOURCE_DIR}   # sources
   ${PROJECT_BINARY_DIR} ) # config.h

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_library( Lib_dsxx ${DRACO_LIBRARY_TYPE} ${sources}  )
if( "${DRACO_LIBRARY_TYPE}" MATCHES "SHARED" )
   set_target_properties( Lib_dsxx 
      PROPERTIES 
      # Provide compile define macro to enable declspec(dllexport) linkage.
      COMPILE_DEFINITIONS BUILDING_DLL 
      # Use custom library naming
      # OUTPUT_NAME ${CMAKE_SHARED_LIBRARY_PREFIX}rtt_ds++${CMAKE_SHARED_LIBRARY_SUFFIX}
      OUTPUT_NAME rtt_ds++
      )
else()
   set_target_properties( Lib_dsxx
      PROPERTIES 
        # Use custom library naming
        # OUTPUT_NAME ${CMAKE_STATIC_LIBRARY_PREFIX}rtt_ds++${CMAKE_STATIC_LIBRARY_SUFFIX}
        OUTPUT_NAME rtt_ds++
      )
endif()

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_dsxx DESTINATION lib )
install( FILES ${headers} DESTINATION include/ds++ )

# ---------------------------------------------------------------------------- #
# Unit tests
# ---------------------------------------------------------------------------- #

if( DRACO_BUILD_TESTS )
   add_subdirectory( test )
endif()   
