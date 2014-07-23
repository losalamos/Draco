#-----------------------------*-cmake-*----------------------------------------#
# file   config/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for diagnostics.
# note   Copyright (C) 2010-2013, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.0.0)
project( diagnostics CXX )

##---------------------------------------------------------------------------##
# Special options for this component
##---------------------------------------------------------------------------##

# DRACO_DIAGNOSICS IS SET IN DS++ (ds++/config.h).

# From autoconf build system: "yes" --> 2, "no" --> 0
set( DRACO_TIMING "0"
   CACHE STRING "Timing detail level (0-2). Default is 0." )
# provide a constrained drop down menu for cmake-gui
set_property( CACHE DRACO_TIMING PROPERTY STRINGS 0 1 2 )

# Are the memory features of getrusage to be activated?
option( USE_PROCMON "Use PROCMON process (cpu/memory) monitoring." OFF )
# If USE_PROCMON is requested but /proc/meminfo is unavailable, abort
if( USE_PROCMON )
   if( NOT EXISTS "/proc/meminfo" )
      message( FATAL_ERROR "USE_PROCMON=ON, but /proc/meminfo unavailable." )
   endif()
endif()
if( NOT ${HAVE_SYS_TYPES} )
   message( FATAL_ERROR "diagnostics/Procmon.hh requires sys/types.h.")
endif()

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
list( REMOVE_ITEM headers
   ${template_implementations} 
   ${PROJECT_SOURCE_DIR}/draco_info_main.cc
   )

# Make the header files available in the IDE.
if( MSVC_IDE OR ${CMAKE_GENERATOR} MATCHES Xcode )
   list( APPEND sources ${headers} )
endif()

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories(
   ${PROJECT_SOURCE_DIR}  # sources
   ${PROJECT_BINARY_DIR}  # config.h
   ${draco_src_dir_SOURCE_DIR} # ds++ and c4 headers
   ${dsxx_BINARY_DIR}     # ds++/config.h
   ${c4_BINARY_DIR}       # c4/config.h
   ${MPI_CXX_INCLUDE_PATH}
   )

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library( 
   TARGET       Lib_diagnostics
   TARGET_DEPS  "Lib_c4"
   LIBRARY_NAME ${PROJECT_NAME}
   SOURCES      "${sources}" )

add_executable( Exe_draco_info ${PROJECT_SOURCE_DIR}/draco_info_main.cc )
target_link_libraries( Exe_draco_info Lib_diagnostics  )
set_target_properties( Exe_draco_info
   PROPERTIES OUTPUT_NAME draco_info 
   PROJECT_LABEL draco_info
   FOLDER        diagnostics )
# Copy required dlls to the executable directory
if (WIN32)
    add_custom_command( TARGET Exe_draco_info POST_BUILD
           COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:Lib_dsxx> $<TARGET_FILE_DIR:Exe_draco_info>
           COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:Lib_c4> $<TARGET_FILE_DIR:Exe_draco_info> )
endif()

# if( QT_FOUND )
#   add_subdirectory( qt )
# endif()

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_diagnostics
   EXPORT draco-targets
   DESTINATION lib )

# ??? Do I need to put draco_info into a bin/ directory for RPATH to
# work correctly after install?
install( TARGETS Exe_draco_info DESTINATION bin )
install( FILES ${headers}  DESTINATION include/diagnostics )

# ---------------------------------------------------------------------------- #
# Unit tests
# ---------------------------------------------------------------------------- #

if( BUILD_TESTING )
   add_subdirectory( test )
endif()   

# ---------------------------------------------------------------------------- #
# Autodoc
# ---------------------------------------------------------------------------- #

process_autodoc_pages()
