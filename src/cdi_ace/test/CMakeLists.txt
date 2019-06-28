#-----------------------------*-cmake-*----------------------------------------#
# file   cdi_ace/test/CMakeLists.txt
# author B. R. Ryan <brryan@lanl.gov>
# date   2019 Jun 19
# brief  Generate build project files for cdi_ace/test.
# note   Copyright (C) 2016-2019, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
project( cdi_ace_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_sources
  ${PROJECT_SOURCE_DIR}/test_Nuclear_Data.cc )
set( input_files
  ${PROJECT_SOURCE_DIR}/1000.14p )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_scalar_tests(
  SOURCES "${test_sources}"
  DEPS "Lib_cdi_ace"
  LABEL "acecheck")

# ---------------------------------------------------------------------------- #
# Copy input files to test directory (binary_dir)
# ---------------------------------------------------------------------------- #

provide_aux_files(
   FILES    "${input_files}" )

#------------------------------------------------------------------------------#
# End cdi_ace/test/CMakeLists.txt
#------------------------------------------------------------------------------#
