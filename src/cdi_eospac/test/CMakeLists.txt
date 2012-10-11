#-----------------------------*-cmake-*----------------------------------------#
# file   cdi_eospac/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for cdi_eospac/test.
# note   Copyright (C) 2010-2012, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id: CMakeLists.txt 6732 2012-09-05 22:28:18Z kellyt $
#------------------------------------------------------------------------------#
project( cdi_eospac_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_sources
   ${PROJECT_SOURCE_DIR}/tEospac.cc
   ${PROJECT_SOURCE_DIR}/tEospacWithCDI.cc
)

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories( 
   ${PROJECT_SOURCE_DIR}      # headers for tests
   ${PROJECT_SOURCE_DIR}/..   # headers for package
   ${PROJECT_BINARY_DIR}/..   # c4/config.h
)

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_scalar_tests(
   SOURCES "${test_sources}" 
   DEPS    "Lib_cdi;${EOSPAC_LIBRARY}" 
   )
