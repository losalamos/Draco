#-----------------------------*-cmake-*----------------------------------------#
# file   c4/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for c4
# note   Copyright (C) 2016-2019, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.9.0)
project( c4 CXX )

#------------------------------------------------------------------------------#
# Quo related setup that will be saved to c4/config.h
#------------------------------------------------------------------------------#
if( TARGET LIBQUO::libquo )
  set( HAVE_LIBQUO 1 )
endif()

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# many c4 and mpi values are set in config/setupMPI.cmake.
# ---------------------------------------------------------------------------- #
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

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #
if( PAPI_LIBRARY )
  set( c4_extra_VENDORS VENDOR_LIST;PAPI;VENDOR_LIBS;${PAPI_LIBRARY} )
endif()

set( target_deps Lib_dsxx )
if ( "${DRACO_C4}" STREQUAL "MPI")
  if( HAVE_LIBQUO )
    list(APPEND target_deps LIBQUO::libquo)
  endif()
  list(APPEND target_deps MPI::MPI_CXX MPI::MPI_C )
  if (TARGET MPI::MPI_Fortran)
    list(APPEND target_deps MPI::MPI_Fortran )
  endif()
endif()

add_component_library(
   TARGET       Lib_c4
   TARGET_DEPS  ${target_deps}
   LIBRARY_NAME c4
   SOURCES      "${sources}"
   HEADERS      "${headers};${template_implementations}"
   ${c4_extra_VENDORS}
   )
target_include_directories( Lib_c4
  PUBLIC $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}> )

# xthi, ythi binaries
if( NOT APPLE )
  add_subdirectory( bin )
endif()

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #
install( TARGETS Lib_c4 EXPORT draco-targets DESTINATION ${DBSCFGDIR}lib )
install( FILES ${headers} DESTINATION ${DBSCFGDIR}include/c4 )
if( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  install(FILES $<TARGET_PDB_FILE:Lib_c4> DESTINATION ${DBSCFGDIR}lib
    OPTIONAL)
endif()

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

#------------------------------------------------------------------------------#
# End c4/CMakeLists.txt
#------------------------------------------------------------------------------#
