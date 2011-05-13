# ---------------------------------------------------------------------------- #
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
# ---------------------------------------------------------------------------- #
# $Id$
# ---------------------------------------------------------------------------- #

project( cdi_gandolf_test CXX )

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

# source_group( cdi_gandolf\\test FILES ${test_lib_sources} ${test_sources} )


