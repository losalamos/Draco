# $Id$

cmake_minimum_required(VERSION 2.6)
project( diagnostics CXX )

# From autoconf build system: "yes" --> 7, "no" --> 0
set( DRACO_DIAGNOSTICS "0"
   CACHE STRING "Diagnostics verbosity level (0-7). Default is 0." )
# From autoconf build system: "yes" --> 2, "no" --> 0
set( DRACO_TIMING "0"
   CACHE STRING "Timing detail level (0-2). Default is 0." )

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #

configure_file( config.h.in ${PROJECT_BINARY_DIR}/diagnostics/config.h )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB template_implementations *.t.hh *.i.hh )
file( GLOB sources *.cc )
file( GLOB explicit_instantiations *_pt.cc )
file( GLOB headers *.hh )
list( APPEND headers ${PROJECT_BINARY_DIR}/diagnostics/config.h )
#list( REMOVE_ITEM headers ${template_implementations} )

# Make the header files available in the IDE.
if( MSVC_IDE )
   list( APPEND sources ${headers} )
endif()

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories(
   ${PROJECT_SOURCE_DIR}  # sources
   ${PROJECT_BINARY_DIR}  # config.h
   ${draco_src_dir_SOURCE_DIR} # ds++ and c4 headers
   ${dsxx_BINARY_DIR}          # ds++/config.h
   ${MPI_INCLUDE_PATH} )

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_library( Lib_diagnostics ${DRACO_LIBRARY_TYPE} ${sources}  )
add_dependencies( Lib_diagnostics Lib_c4 Lib_dsxx )
if( "${DRACO_LIBRARY_TYPE}" MATCHES "SHARED" )
   set_target_properties( Lib_diagnostics 
      PROPERTIES 
      # Provide compile define macro to enable declspec(dllexport) linkage.
      COMPILE_DEFINITIONS BUILDING_DLL 
      # Use custom library naming
      # OUTPUT_NAME ${CMAKE_SHARED_LIBRARY_PREFIX}rtt_ds++${CMAKE_SHARED_LIBRARY_SUFFIX}
      OUTPUT_NAME rtt_diagnostics
      )
else()
   set_target_properties( Lib_diagnostics
      PROPERTIES 
        # Use custom library naming
        # OUTPUT_NAME ${CMAKE_STATIC_LIBRARY_PREFIX}rtt_ds++${CMAKE_STATIC_LIBRARY_SUFFIX}
        OUTPUT_NAME rtt_diagnostics
      )
endif()

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_diagnostics DESTINATION lib )
install( FILES ${headers}  DESTINATION include/diagnostics )

# ---------------------------------------------------------------------------- #
# Unit tests
# ---------------------------------------------------------------------------- #

if( DRACO_BUILD_TESTS )
   add_subdirectory( test )
endif()   
