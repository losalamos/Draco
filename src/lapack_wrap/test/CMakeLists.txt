# ---------------------------------------------------------------------------- #
# Draco - Support libraries for transport solvers.
# Copyright (C) 2016-2017 Los Alamos National Security, LLC.
# All rights reserved.
# ---------------------------------------------------------------------------- #
# $Id$
# ---------------------------------------------------------------------------- #

project( lapack_wrap_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB test_sources *.cc )

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories(
   ${PROJECT_SOURCE_DIR}      # headers for tests
   ${PROJECT_SOURCE_DIR}/..   # headers for package
)

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

set( test_dependencies
   Lib_dsxx
   lapack
   blas
   )

add_scalar_tests(
   SOURCES "${test_sources}"
   DEPS    "${test_dependencies}" )
