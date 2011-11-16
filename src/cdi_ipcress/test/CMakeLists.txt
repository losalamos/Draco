# ---------------------------------------------------------------------------- #
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
# ---------------------------------------------------------------------------- #
# $Id$
# ---------------------------------------------------------------------------- #

project( cdi_ipcress_test CXX )

# ---------------------------------------------------------------------------- #
# Configuration files
# ---------------------------------------------------------------------------- #
set( IPCRESS_INTERPRETER_BIN_DIR ${cdi_ipcress_BINARY_DIR}/bin )
set( IPCRESS_INTERPRETER_SOURCE_DIR ${PROJECT_SOURCE_DIR} )
set( IPCRESS_INTERPRETER_BUILD_DIR ${PROJECT_BINARY_DIR} )

configure_file( config.h.in ${PROJECT_BINARY_DIR}/config.h )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_lib_sources 
   ${PROJECT_SOURCE_DIR}/cdi_ipcress_test.cc 
   ${PROJECT_SOURCE_DIR}/cdi_ipcress_test.hh
)
# file( GLOB test_sources *.cc )
# list( REMOVE_ITEM test_sources ${test_lib_sources} )
set( test_sources 
   tIpcressFile.cc
   tIpcressOpacity.cc
   tIpcressWithCDI.cc
   tIpcress_Interpreter.cc
)
file( GLOB ipcress_data *.ipcress )

##---------------------------------------------------------------------------##
## Copy data files
##---------------------------------------------------------------------------##

foreach( file ${ipcress_data} )
   get_filename_component( shortname ${file} NAME )
   configure_file( ${file} 
      ${PROJECT_BINARY_DIR}/${shortname} COPYONLY )
endforeach()

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories( 
   ${PROJECT_SOURCE_DIR}      # headers for tests
   ${PROJECT_BINARY_DIR}      # test/config.h
   ${PROJECT_SOURCE_DIR}/..   # headers for package
)

# ---------------------------------------------------------------------------- #
# Build library for test directory
# ---------------------------------------------------------------------------- #

add_component_library( Lib_cdi_ipcress_test ${PROJECT_NAME} ${test_lib_sources} )
add_dependencies(      Lib_cdi_ipcress_test Lib_cdi_ipcress )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_scalar_tests(
   SOURCES "${test_sources}"
   DEPS    "Lib_cdi;Lib_dsxx" )



