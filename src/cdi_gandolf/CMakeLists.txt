# $Id$

cmake_minimum_required(VERSION 2.6)
project( cdi_gandolf CXX )

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #

# 2010-12-3 KT: No longer supporting stubbed out cdi_ganodolf
# if( NOT GANDOLF_FOUND )
#    set( rtt_cdi_gandolf_stub 1 )
# endif()

configure_file( config.h.in ${PROJECT_BINARY_DIR}/cdi_gandolf/config.h )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB sources *.cc )
file( GLOB explicit_instantiations *_pt.cc )
file( GLOB headers *.hh )
list( APPEND headers ${PROJECT_BINARY_DIR}/cdi_gandolf/config.h )

# Make the header files available in the IDE.
if( MSVC_IDE OR ${CMAKE_GENERATOR} MATCHES Xcode )
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

# setup binary directory
file( MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/bin )

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library( Lib_cdi_gandolf ${PROJECT_NAME} "${sources}" )
add_dependencies(      Lib_cdi_gandolf Lib_cdi Lib_dsxx )


add_executable( Exe_Gandolf_Interpreter ${sources} )
target_link_libraries( Exe_Gandolf_Interpreter 
                       Lib_cdi_gandolf Lib_cdi Lib_dsxx 
                       ${GANDOLF_LIBRARY} ${GANDOLF_GFORTRAN_LIBRARY})
set_target_properties( Exe_Gandolf_Interpreter 
                       PROPERTIES OUTPUT_NAME bin/Gandolf_Interpreter )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_cdi_gandolf DESTINATION lib )
install( FILES ${headers} DESTINATION include/cdi_gandolf )
install( TARGETS Exe_Gandolf_Interpreter DESTINATION bin )

# ---------------------------------------------------------------------------- #
# Unit tests
# ---------------------------------------------------------------------------- #

if( BUILD_TESTING AND EXISTS ${PROJECT_SOURCE_DIR}/test/CMakeLists.txt )
   add_subdirectory( test )
endif()   

