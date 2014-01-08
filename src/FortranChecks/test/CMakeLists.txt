# ---------------------------------------------------------------------------- #
# file   FortranChecks/test/CMakeLists.txt
# author Kelly Thompson
# brief  Instructions for building a C++ main that links against an
#        F90 library and uses ISO_C_BINDING.
# note   Copyright (C) 2010-2013 Los Alamos National Security, LLC.
#        All rights reserved
# ---------------------------------------------------------------------------- #
# $Id$
# ---------------------------------------------------------------------------- #

project( FortranChecks_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB test_sources *.cc )

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories( 
   ${PROJECT_SOURCE_DIR}      # headers for tests
   ${dsxx_SOURCE_DIR}/..      # headers for ds++
   ${dsxx_BINARY_DIR}   )     # ds++/config.h

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

# Create the test.
add_scalar_tests( 
   SOURCES "${test_sources}" 
   DEPS    "Lib_dsxx;Lib_FortranChecks" )

# ---------------------------------------------------------------------------- #
# end of FortranChecks/test/CMakeLists.txt
# ---------------------------------------------------------------------------- #
