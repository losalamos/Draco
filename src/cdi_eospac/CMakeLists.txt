#-----------------------------*-cmake-*----------------------------------------#
# file   src/cdi_eospac/CMakeLists.txt
# brief  Instructions for building cdi_eospac Makefiles.
# note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.0.0)
project( cdi_eospac CXX )

##---------------------------------------------------------------------------##
# Requires EOSPAC
##---------------------------------------------------------------------------##
find_package( EOSPAC QUIET )
set_package_properties( EOSPAC PROPERTIES
   DESCRIPTION "Access SESAME thermodynamic and transport data."
   TYPE OPTIONAL
   PURPOSE "Required for bulding the cdi_eospac component."
   )

if( EOSPAC_FOUND )
# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB sources *.cc )
file( GLOB headers *.hh )
list( REMOVE_ITEM sources ${PROJECT_SOURCE_DIR}/QueryEospac.cc )

# Make the header files available in the IDE.
if( MSVC_IDE OR ${CMAKE_GENERATOR} MATCHES Xcode )
   list( APPEND sources ${headers} )
endif()

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories(
   ${PROJECT_SOURCE_DIR}      # component headers
   ${PROJECT_BINARY_DIR}      # cdi_eospac/config.h
   ${PROJECT_SOURCE_DIR}/..   # ds++ and cdi header files
   ${dsxx_BINARY_DIR}         # ds++/config.h
   ${ode_BINARY_DIR}          # ode/config.h
   ${EOSPAC_INCLUDE_DIR}
)

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library(
   TARGET       Lib_cdi_eospac
   TARGET_DEPS  Lib_cdi
   LIBRARY_NAME cdi_eospac
   SOURCES      "${sources}"
   VENDOR_LIST "EOSPAC"
   VENDOR_LIBS "${EOSPAC_LIBRARY}"
   VENDOR_INCLUDE_DIRS "${EOSPAC_INCLUDE_DIR}" )

add_component_executable(
  TARGET      Exe_QueryEospac
  TARGET_DEPS Lib_cdi_eospac
  SOURCES     ${PROJECT_SOURCE_DIR}/QueryEospac.cc
  PREFIX       Draco
  )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_cdi_eospac  EXPORT draco-targets DESTINATION ${DBSCFGDIR}lib )
install( TARGETS Exe_QueryEospac EXPORT draco-targets DESTINATION ${DBSCFGDIR}bin )
install( FILES ${headers} DESTINATION ${DBSCFGDIR}include/cdi_eospac )

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

endif() # EOSPAC_FOUND

# ---------------------------------------------------------------------------- #
# End cdi_eospac/CMakeLists.txt
# ---------------------------------------------------------------------------- #