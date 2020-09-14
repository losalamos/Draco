#--------------------------------------------*-cmake-*---------------------------------------------#
# file   c4/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for c4
# note   Copyright (C) 2016-2020, Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.17.0)
project( c4 CXX )

#--------------------------------------------------------------------------------------------------#
# Quo related setup that will be saved to c4/config.h
#--------------------------------------------------------------------------------------------------#
if( TARGET LIBQUO::libquo AND DRACO_C4 STREQUAL MPI )
  set( HAVE_LIBQUO 1 )
endif()

#--------------------------------------------------------------------------------------------------#
# Generate config.h (only occurs when cmake is run) many c4 and MPI values are set in
# config/setupMPI.cmake.
#--------------------------------------------------------------------------------------------------#
configure_file( config.h.in ${PROJECT_BINARY_DIR}/c4/config.h )

#--------------------------------------------------------------------------------------------------#
# Source files
#--------------------------------------------------------------------------------------------------#
file( GLOB sources *.cc )
file( GLOB headers *.hh *.h)
list( APPEND headers ${PROJECT_BINARY_DIR}/c4/config.h )

#--------------------------------------------------------------------------------------------------#
# Build package library
#--------------------------------------------------------------------------------------------------#
if( PAPI_LIBRARY )
  set( c4_extra_VENDORS VENDOR_LIST;PAPI;VENDOR_LIBS;${PAPI_LIBRARY} )
endif()

set( target_deps Lib_dsxx )
if ( "${DRACO_C4}" STREQUAL "MPI")
  if( HAVE_LIBQUO )
    list(APPEND target_deps LIBQUO::libquo)
  endif()
  list(APPEND target_deps MPI::MPI_CXX MPI::MPI_C )
endif()

add_component_library(
   TARGET       Lib_c4
   TARGET_DEPS  ${target_deps}
   LIBRARY_NAME c4
   SOURCES      "${sources}"
   HEADERS      "${headers}"
   ${c4_extra_VENDORS}
   )
target_include_directories( Lib_c4
  PUBLIC $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}> )

#--------------------------------------------------------------------------------------------------#
# Fortran-based MPI wrapper
#--------------------------------------------------------------------------------------------------#

if( HAVE_Fortran )

  set( deplist Lib_c4 Lib_dsxx )

  # For Unix, this will simply run add_subdirectory(fc4). For Visual Studio or XCode, an
  # ExternalProject will be created that uses an alternate gfortran via Makefiles to generate the
  # library external to the main project.
  include(CMakeAddFortranSubdirectory)
  init_build_system_state() # returns ${build_system_state} used below.

  cmake_add_fortran_subdirectory(
    fc4                       # directory name
    PROJECT c4_fc4            # project name in top level CMakeLists.txt
                              # Creates target named 'c4_fc4_build'
    ARCHIVE_DIR  fc4          # .lib location relative to root binary tree
    RUNTIME_DIR  fc4          # .dll location relative to root binary tree
    LIBRARIES    rtt_fc4      # new library will have this name. --> librtt_fc4.dll.
    TARGET_NAMES Lib_c4_fc4
    DEPENDS      "${deplist}"
    # NO_EXTERNAL_INSTALL - we want to install this target!
    #                     - if not set, the cafs function will set CAFS_EXPORT_TARGET_PROPERTIES and
    #                       CAFS_EXPORT_DEFINE_IMPORT_PREFIX
    CMAKE_COMMAND_LINE "${build_system_state}"
    # VERBOSE
  )
  string(APPEND Draco_EXPORT_TARGET_PROPERTIES
    "${CAFS_EXPORT_DEFINE_IMPORT_PREFIX}"
    "${CAFS_EXPORT_TARGET_PROPERTIES}")
  set(Draco_EXPORT_TARGET_PROPERTIES "${Draco_EXPORT_TARGET_PROPERTIES}"
    PARENT_SCOPE)

  cmake_add_fortran_subdirectory(
    ftest                       # directory name
    PROJECT c4_ftest            # project name in top level CMakeLists.txt
                                # Creates target named 'quadrature_ftest_build'
    ARCHIVE_DIR  ftest          # .lib location relative to root binary tree
    RUNTIME_DIR  ftest          # .dll location relative to root binary tree
    LIBRARIES    rtt_c4_ftest   # new library will have this name. --> librtt_c4_ftest.dll.
    TARGET_NAMES Lib_c4_ftest
    DEPENDS      "Lib_c4_fc4;${deplist}"
    NO_EXTERNAL_INSTALL
    CMAKE_COMMAND_LINE "${build_system_state}"
    # VERBOSE
    )

endif()

# xthi, ythi binaries
add_subdirectory( bin )

#--------------------------------------------------------------------------------------------------#
# Installation instructions
#--------------------------------------------------------------------------------------------------#
install( TARGETS Lib_c4 EXPORT draco-targets DESTINATION ${DBSCFGDIR}lib )
install( FILES ${headers} DESTINATION ${DBSCFGDIR}include/c4 )
if( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  install(FILES $<TARGET_PDB_FILE:Lib_c4> DESTINATION ${DBSCFGDIR}lib
    OPTIONAL)
endif()

#--------------------------------------------------------------------------------------------------#
# Unit tests
#--------------------------------------------------------------------------------------------------#
if( BUILD_TESTING )
  add_subdirectory( test )
endif()

#--------------------------------------------------------------------------------------------------#
# Autodoc
#--------------------------------------------------------------------------------------------------#
process_autodoc_pages()

#--------------------------------------------------------------------------------------------------#
# End c4/CMakeLists.txt
#--------------------------------------------------------------------------------------------------#
