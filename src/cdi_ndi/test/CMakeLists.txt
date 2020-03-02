#-----------------------------*-cmake-*----------------------------------------#
# file   cdi_ndi/test/CMakeLists.txt
# author Ben. R. Ryan <brryan@lanl.gov>
# date   2020 Feb 4
# brief  Generate build project files for cdi_ndi/test.
# note   Copyright (C) 2020, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
project( cdi_ndi_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_sources
  tstNDI_TN.cc )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_scalar_tests(
  SOURCES "${test_sources}"
  DEPS    "Lib_cdi_ndi" )

#------------------------------------------------------------------------------#
# End cdi_ndi/test/CMakeLists.txt
#------------------------------------------------------------------------------#
