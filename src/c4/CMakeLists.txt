# $Id$

cmake_minimum_required(VERSION 2.6)
project( c4 CXX )

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #

if( "${DRACO_C4}" MATCHES "MPI" )
  set( C4_MPI 1 CACHE STRING "Are we building an MPI aware version?"  )
else()
  set( C4_SCALAR 1 CACHE STRING "Are we building a scalar-only version (no mpi)?" )
endif()
if( UNIX )
  set( c4_isLinux 1 )
elseif( WIN32 )
  set( c4_isWin 1 )
elseif( OF1 )
  set( c4_isOSF1 1 )
elseif( APPLE )
  set( c4_isDarwin 1)
else()
  set( c4_isAIX 1 )
endif()
configure_file( config.h.in ${PROJECT_BINARY_DIR}/c4/config.h )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB template_implementations *.t.hh *.i.hh )
file( GLOB sources *.cc )
file( GLOB explicit_instantiations *_pt.cc )
file( GLOB headers *.hh *.h)
list( APPEND headers ${PROJECT_BINARY_DIR}/c4/config.h )
list( REMOVE_ITEM headers ${template_implementations} )

# Make the header files available in the IDE.
if( MSVC_IDE )
   list( APPEND sources ${headers} )
endif()

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories( ${PROJECT_SOURCE_DIR}       # sources
                     ${PROJECT_BINARY_DIR}       # c4/config.h
                     ${draco_src_dir_SOURCE_DIR} # ds++ header files
                     ${dsxx_BINARY_DIR}          # ds++/config.h
                     ${MPI_INCLUDE_PATH}
)

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_library( Lib_c4 ${DRACO_LIBRARY_TYPE} ${sources}  )
#add_dependencies( Lib_c4 Lib_dsxx )
if( "${DRACO_LIBRARY_TYPE}" MATCHES "SHARED" )
  set_target_properties( Lib_c4 
    PROPERTIES 
      # Provide compile define macro to enable declspec(dllexport) linkage.
      COMPILE_DEFINITIONS BUILDING_DLL 
      # Use custom library naming
      OUTPUT_NAME rtt_c4
  )
else()
  set_target_properties( Lib_c4 
    PROPERTIES 
      # Use custom library naming
      OUTPUT_NAME rtt_c4
  )
endif()

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_c4 
         DESTINATION lib 
         EXPORT draco_c4 )
install( FILES ${headers} DESTINATION include/c4 )
install( EXPORT draco_c4
         DESTINATION . )

# ---------------------------------------------------------------------------- #
# Unit tests
# ---------------------------------------------------------------------------- #

if( DRACO_BUILD_TESTS )
 add_subdirectory( test )
endif()   
