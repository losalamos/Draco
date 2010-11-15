# $Id$

cmake_minimum_required(VERSION 2.6)
project( RTT_Format_Reader CXX )

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
                     ${MPI_INCLUDE_PATH}
)

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_library( Lib_RTT_Format_Reader ${DRACO_LIBRARY_TYPE} ${sources} )
add_dependencies( Lib_RTT_Format_Reader
   Lib_meshReaders
   Lib_mesh_element 
   Lib_dsxx )
if( "${DRACO_LIBRARY_TYPE}" MATCHES "SHARED" )
  set_target_properties( Lib_RTT_Format_Reader 
    PROPERTIES 
      # Provide compile define macro to enable declspec(dllexport) linkage.
      COMPILE_DEFINITIONS BUILDING_DLL 
      # Use custom library naming
      OUTPUT_NAME rtt_RTT_Format_Reader
  )
else()
  set_target_properties( Lib_RTT_Format_Reader 
    PROPERTIES 
      # Use custom library naming
      OUTPUT_NAME rtt_RTT_Format_Reader
  )
endif()

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_RTT_Format_Reader DESTINATION lib )
install( FILES ${headers} DESTINATION include/RTT_Format_Reader )

# ---------------------------------------------------------------------------- #
# Unit tests
# ---------------------------------------------------------------------------- #

if( DRACO_BUILD_TESTS )
 add_subdirectory( test )
endif()   
