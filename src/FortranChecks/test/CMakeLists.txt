# ---------------------------------------------------------------------------- #
# file   FortranChecks/test/CMakeLists.txt
# author Kelly Thompson
# brief  Instructions for building a C++ main that links against an
#        F90 library and uses ISO_C_BINDING.
# note   Copyright (C) 2010-2014 Los Alamos National Security, LLC.
#        All rights reserved
# ---------------------------------------------------------------------------- #
# $Id$
#------------------------------------------------------------------------------#
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

# When building on SQ with XLC/XLF90, ensure that libxlomp_ser is not
# on the link line.  This macro is defined in config/unix-xlf.cmake.
if( ${CMAKE_Fortran_COMPILER_FLAVOR} MATCHES "XL" )
   remove_lib_from_link( "xlomp_ser" )
endif()

# Create the test.
set( test_deps "Lib_dsxx;Lib_FC_f90sub" )
add_scalar_tests( 
   SOURCES "${test_sources}" 
   DEPS    "${test_deps}" )

# ---------------------------------------------------------------------------- #
# end of FortranChecks/test/CMakeLists.txt
# ---------------------------------------------------------------------------- #

