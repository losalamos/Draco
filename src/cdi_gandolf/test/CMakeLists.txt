# ---------------------------------------------------------------------------- #
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
# ---------------------------------------------------------------------------- #
# $Id$
# ---------------------------------------------------------------------------- #

project( cdi_gandolf_test CXX )

# ---------------------------------------------------------------------------- #
# Configuration files
# ---------------------------------------------------------------------------- #
set( GANDOLF_INTERPRETER_BIN_DIR ${cdi_gandolf_BINARY_DIR}/bin )
set( GANDOLF_INTERPRETER_SOURCE_DIR ${PROJECT_SOURCE_DIR} )
set( GANDOLF_INTERPRETER_BUILD_DIR ${PROJECT_BINARY_DIR} )

configure_file( config.h.in ${PROJECT_BINARY_DIR}/config.h )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_lib_sources 
   ${PROJECT_SOURCE_DIR}/cdi_gandolf_test.cc 
   ${PROJECT_SOURCE_DIR}/cdi_gandolf_test.hh
)
file( GLOB test_sources *.cc )
list( REMOVE_ITEM test_sources ${test_lib_sources} )
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

add_component_library( Lib_cdi_gandolf_test ${PROJECT_NAME} ${test_lib_sources} )
add_dependencies(      Lib_cdi_gandolf_test Lib_cdi_gandolf )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_scalar_tests(
   SOURCES "${test_sources}"
   DEPS    "Lib_cdi;Lib_dsxx;${GANDOLF_LIBRARIES}" )



