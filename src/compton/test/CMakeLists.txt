#-----------------------------*-cmake-*----------------------------------------#
# file   compton/test/CMakeLists.txt
# author Kendra Keady
# date   2017 February 28
# brief  Generate build project files for compton/test.
# note   Copyright (C) 2017-2020, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

project( compton_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_sources
  ${PROJECT_SOURCE_DIR}/tCompton.cc )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

# If the compton library was not found, then there is no Lib_compton library to
# link against (and thus compton/config.h isn't found). However, we can still
# build the test against the stub headers...
if( NOT TARGET Lib_compton )

  add_parallel_tests(
    SOURCES "${test_sources}"
    PE_LIST "1"
    DEPS    "Lib_c4"
    LABEL   "nomemcheck"
    MPI_PLUS_OMP )

  # generated include directive files (config.h)
  target_include_directories( Ut_compton_tCompton_exe
    PUBLIC $<BUILD_INTERFACE:${compton_BINARY_DIR}> )

else()

  add_parallel_tests(
    SOURCES "${test_sources}"
    PE_LIST "1"
    DEPS    "Lib_compton;Lib_c4"
    LABEL   "nomemcheck"
    MPI_PLUS_OMP )

endif()

# ---------------------------------------------------------------------------- #
# end of compton/test/ CMakeLists.txt
# ---------------------------------------------------------------------------- #
