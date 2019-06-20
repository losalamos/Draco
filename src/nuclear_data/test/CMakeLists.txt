#-----------------------------*-cmake-*----------------------------------------#
# file   nuclear_data/test/CMakeLists.txt
# author B. R. Ryan <brryan@lanl.gov>
# date   2019 Jun 19
# brief  Generate build project files for nuclear_data/test.
# note   Copyright (C) 2016-2019, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
project( nuclear_data_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB test_sources *.cc )
#file( GLOB mesh_input_files *.mesh.in )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

#add_scalar_tests(
#   SOURCES "${test_sources}"
#   DEPS    "Lib_nuclear_data" )

# ---------------------------------------------------------------------------- #
# Copy input files to test directory (binary_dir)
# ---------------------------------------------------------------------------- #

#provide_aux_files(
#   FILES    "${mesh_input_files}"
#   SRC_EXT  ".mesh.in"
#   DEST_EXT ".mesh" )

#------------------------------------------------------------------------------#
# End nuclear_data/test/CMakeLists.txt
#------------------------------------------------------------------------------#
