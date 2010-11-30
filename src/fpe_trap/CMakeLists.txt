# $Id$

cmake_minimum_required(VERSION 2.6)
project( fpe_trap CXX )

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #

if( ${CMAKE_SYSTEM} MATCHES "OSF" )
   set( FPETRAP_OSF_ALPHA 1 )
elseif( ${CMAKE_SYSTEM} MATCHES "Linux" )
   set( FPETRAP_LINUX_X86 1 )
elseif( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
   if( ${_CMAKE_OSX_MACHINE} MATCHES "ppc" )
      set( FPETRAP_DARWIN_PPC 1 )
   else()
      set( FPETRAP_DARWIN_INTEL 1 )
   endif()
endif()

configure_file( config.h.in ${PROJECT_BINARY_DIR}/fpe_trap/config.h )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

# file( GLOB template_implementations *.t.hh *.i.hh )
file( GLOB sources *.cc )
file( GLOB explicit_instantiations *_pt.cc )
file( GLOB headers *.hh )
list( APPEND headers ${PROJECT_BINARY_DIR}/fpe_trap/config.h )
#list( REMOVE_ITEM headers ${template_implementations} )

# Make the header files available in the IDE.
if( MSVC_IDE )
   list( APPEND sources ${headers} )
endif()

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories( 
   ${PROJECT_SOURCE_DIR}      # sources
   ${PROJECT_BINARY_DIR}      # config.h
   ${PROJECT_SOURCE_DIR}/..   # ds++ header files
   ${dsxx_BINARY_DIR}         # ds++/config.h
) 

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_library( Lib_fpe_trap ${DRACO_LIBRARY_TYPE} ${sources}  )
add_dependencies( Lib_fpe_trap Lib_cdi Lib_dsxx )
if( "${DRACO_LIBRARY_TYPE}" MATCHES "SHARED" )
   set_target_properties( Lib_fpe_trap 
      PROPERTIES 
      # Provide compile define macro to enable declspec(dllexport) linkage.
      COMPILE_DEFINITIONS BUILDING_DLL 
      # Use custom library naming
      OUTPUT_NAME rtt_fpe_trap
      )
else()
   set_target_properties( Lib_fpe_trap 
      PROPERTIES 
      # Use custom library naming
      OUTPUT_NAME rtt_fpe_trap
      )
endif()

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_fpe_trap DESTINATION lib )
install( FILES ${headers} DESTINATION include/fpe_trap )

# ---------------------------------------------------------------------------- #
# Unit tests
# ---------------------------------------------------------------------------- #

if( DRACO_BUILD_TESTS )
   add_subdirectory( test )
endif()   
