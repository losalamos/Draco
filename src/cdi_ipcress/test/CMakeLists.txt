#-----------------------------*-cmake-*----------------------------------------#
# file   cdi_ipcress/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for cdi_ipcress/test.
# note   Copyright (C) 2010-2015, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#
project( cdi_ipcress_test CXX )

# ---------------------------------------------------------------------------- #
# Configuration files
# ---------------------------------------------------------------------------- #

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
set( test_sources
   tIpcressFile.cc
   tIpcressOpacity.cc
   tIpcressWithCDI.cc
   tIpcress_Interpreter.cc
#   ReadOdfIpcressFile.cc # needs argument '--unittest'
   tIpcressOdfmgOpacity.cc
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

add_component_library(
   TARGET       Lib_cdi_ipcress_test
   TARGET_DEPS  Lib_cdi_ipcress
   LIBRARY_NAME ${PROJECT_NAME}
   SOURCES      "${test_lib_sources}" )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_scalar_tests(
   SOURCES "${test_sources}"
   DEPS    "Lib_cdi_ipcress_test" )

add_scalar_tests(
   SOURCES "ReadOdfIpcressFile.cc"
   DEPS    "Lib_cdi_ipcress_test"
   TEST_ARGS "--unittest" )