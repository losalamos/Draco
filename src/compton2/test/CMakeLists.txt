#-----------------------------*-cmake-*----------------------------------------#
# file   compton2/test/CMakeLists.txt
# author Andrew Till
# date   11 May 2020
# brief  Generate build project files for compton2/test.
# note   Copyright (C) 2017-2020, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

project( compton2_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_sources
  ${PROJECT_SOURCE_DIR}/tCompton2.cc )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

if( NOT TARGET Lib_compton2 )

  add_parallel_tests(
    SOURCES "${test_sources}"
    PE_LIST "1"
    DEPS    "Lib_c4"
    LABEL   "nomemcheck"
    MPI_PLUS_OMP )

  # generated include directive files (config.h)
  target_include_directories( Ut_compton_tCompton2_exe
    PUBLIC $<BUILD_INTERFACE:${compton_BINARY_DIR}> )

else()

  add_parallel_tests(
    SOURCES "${test_sources}"
    PE_LIST "1"
    DEPS    "Lib_compton2;Lib_c4"
    LABEL   "nomemcheck"
    MPI_PLUS_OMP )

endif()

# ---------------------------------------------------------------------------- #
# end of compton2/test/CMakeLists.txt
# ---------------------------------------------------------------------------- #
