# $Id$

cmake_minimum_required(VERSION 2.6)
project( quadrature CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

#file( GLOB template_implementations *.t.hh *.i.hh )
file( GLOB sources *.cc )
#file( GLOB explicit_instantiations *_pt.cc )
file( GLOB headers *.hh )
#list( REMOVE_ITEM headers ${template_implementations} )

# Make the header files available in the IDE.
if( MSVC_IDE )
   list( APPEND sources ${headers} )
endif()

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories( ${PROJECT_SOURCE_DIR}       # sources
                     ${draco_src_dir_SOURCE_DIR} # ds++ header files
                     ${dsxx_BINARY_DIR}          # ds++/config.h
                     ${GSL_INCLUDE_DIRS}
                     ${MPI_INCLUDE_PATH}
)

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_library( Lib_quadrature ${DRACO_LIBRARY_TYPE} ${sources} )
add_dependencies( Lib_quadrature
   Lib_units
   Lib_special_functions
   Lib_ode )
if( "${DRACO_LIBRARY_TYPE}" MATCHES "SHARED" )
  set_target_properties( Lib_quadrature 
    PROPERTIES 
      # Provide compile define macro to enable declspec(dllexport) linkage.
      COMPILE_DEFINITIONS BUILDING_DLL 
      # Use custom library naming
      OUTPUT_NAME rtt_quadrature
  )
else()
  set_target_properties( Lib_quadrature 
    PROPERTIES 
      # Use custom library naming
      OUTPUT_NAME rtt_quadrature
  )
endif()

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_quadrature DESTINATION lib )
install( FILES ${headers} DESTINATION include/quadrature )

# ---------------------------------------------------------------------------- #
# Unit tests
# ---------------------------------------------------------------------------- #

if( DRACO_BUILD_TESTS )
 add_subdirectory( test )
endif()   
