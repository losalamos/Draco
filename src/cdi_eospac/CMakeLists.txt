#-----------------------------*-cmake-*----------------------------------------#
# file   src/cdi_eospac/CMakeLists.txt
# brief  Instructions for building cdi_eospac Makefiles.
# note   Copyright (C) 2010-2012 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id: CMakeLists.txt 6732 2012-09-05 22:28:18Z kellyt $
#------------------------------------------------------------------------------#
cmake_minimum_required( VERSION 2.8 )
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
add_executable( Exe_QueryEospac QueryEospac.cc )
target_link_libraries( Exe_QueryEospac Lib_cdi_eospac )
set_target_properties( Exe_QueryEospac
   PROPERTIES OUTPUT_NAME QueryEospac )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_cdi_eospac EXPORT draco-targets DESTINATION lib )
install( FILES ${headers} DESTINATION include/cdi_eospac )
install( TARGETS Exe_QueryEospac DESTINATION bin )

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