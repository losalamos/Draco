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
# Source files
# ---------------------------------------------------------------------------- #

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
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_scalar_tests(
   SOURCES "${test_sources}"
   DEPS    "Lib_cdi_ipcress" )

add_scalar_tests(
   SOURCES "ReadOdfIpcressFile.cc"
   DEPS    "Lib_cdi_ipcress"
   TEST_ARGS "--unittest" )
   
# tIpcress_Interpreter needs to know where the binary Ipcress_Interpreter is located.  
# Since the build system knows this location and it may be difficult to find from the
# unit test, simply provide the location via CPP define.
set_target_properties( Ut_cdi_ipcress_tIpcress_Interpreter_exe
   PROPERTIES COMPILE_DEFINITIONS "IILOC=\"$<TARGET_FILE:Exe_Ipcress_Interpreter>\"" )
