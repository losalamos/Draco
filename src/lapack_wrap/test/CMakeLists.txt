# ---------------------------------------------------------------------------- #
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
# ---------------------------------------------------------------------------- #
# $Id$
# ---------------------------------------------------------------------------- #

project( lapack_wrap_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_lib_sources 
   ${PROJECT_SOURCE_DIR}/lapack_wrap_test.cc 
   ${PROJECT_SOURCE_DIR}/lapack_wrap_test.hh
)
file( GLOB test_sources *.cc )
list( REMOVE_ITEM test_sources ${test_lib_sources} )

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

add_component_library( Lib_lapack_wrap_test ${PROJECT_NAME} 
   "${test_lib_sources}" )
add_dependencies(      Lib_lapack_wrap_test Lib_lapack_wrap )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

set( test_dependencies
   Lib_dsxx
   ${LAPACK_LIBRARIES}
   ${INTEL_LIBM}
   )

# [KT:2011-02-17]  I can't find where this is set
#    ${FORTRAN_LIBS}
#    ${PGI_EXTRA_F90_LIBS}

add_scalar_tests(
   SOURCES "${test_sources}"
   DEPS    "${test_dependencies}" )

