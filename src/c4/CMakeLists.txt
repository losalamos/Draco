#-----------------------------*-cmake-*----------------------------------------#
# file   c4/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for c4
# note   Copyright (C) 2010-2012, Los Alamos National Security
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 2.6)
project( c4 CXX )

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #

if( "${DRACO_C4}" MATCHES "MPI" )
  set( C4_SCALAR 0 )
  set( C4_MPI 1 )
else()
  set( C4_SCALAR 1 )
  set( C4_MPI 0 )
endif()
set( C4_SCALAR ${C4_SCALAR} CACHE STRING 
   "Are we building a scalar-only version (no mpi in c4/config.h)?"
   FORCE )
set( C4_MPI ${C4_MPI} CACHE STRING 
   "Are we building an MPI aware version? (c4/config.h)" FORCE )
mark_as_advanced( C4_MPI C4_SCALAR )

if( UNIX )
   set( C4_MPICMD "mpirun ${MPIEXEC_POSTFLAGS_STRING} -np " )
elseif( WIN32 )
   set( C4_MPICMD "mpiexec -np " )
elseif( OSF1 )
   set( C4_MPICMD "prun -n " )
elseif( APPLE )
   set( C4_MPICMD "mpirun -np " )
else() # AIX
   set( C4_MPICMD "poe -procs " )
endif()
if( "${SITE}" MATCHES "c[it]" )
   set( C4_MPICMD "aprun -n " )
endif()

configure_file( config.h.in ${PROJECT_BINARY_DIR}/c4/config.h )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB template_implementations *.t.hh )
file( GLOB sources *.cc )
file( GLOB explicit_instantiations *_pt.cc )
file( GLOB headers *.hh *.h)
list( APPEND headers ${PROJECT_BINARY_DIR}/c4/config.h )
list( REMOVE_ITEM headers ${template_implementations} )

# Make the header files available in the IDE.
if( MSVC_IDE OR ${CMAKE_GENERATOR} MATCHES Xcode )
   list( APPEND sources ${headers} )
endif()

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories( ${PROJECT_SOURCE_DIR}       # sources
                     ${PROJECT_BINARY_DIR}       # c4/config.h
                     ${draco_src_dir_SOURCE_DIR} # ds++ header files
                     ${dsxx_BINARY_DIR}          # ds++/config.h
                     )

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

set( c4_VENDOR_LIBS "${MPI_LIBRARIES}")
if( PAPI_LIBRARY )
   list( APPEND c4_VENDOR_LIBS ${PAPI_LIBRARY} )
endif()

add_component_library( 
   TARGET       Lib_c4
   TARGET_DEPS  Lib_dsxx
   LIBRARY_NAME c4 
   SOURCES      "${sources}"
   VENDOR_LIST  "MPI" 
   VENDOR_LIBS  "${c4_VENDOR_LIBS}"
   VENDOR_INCLUDE_DIRS "${MPI_CXX_INCLUDE_PATH}" 
   )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_c4 EXPORT draco-targets DESTINATION lib )
install( FILES ${headers} DESTINATION include/c4 )

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
