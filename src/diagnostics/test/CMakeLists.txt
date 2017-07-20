#-----------------------------*-cmake-*----------------------------------------#
# file   diagnostics/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for diagnostics/test.
# note   Copyright (C) 2016-2017, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
project( diagnostics_test CXX )

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
   ${PROJECT_BINARY_DIR}/..   # config.h
   )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #
add_scalar_tests(
   SOURCES "${test_sources}"
   DEPS    Lib_diagnostics )

include( ApplicationUnitTest )

add_app_unit_test(
  DRIVER ${CMAKE_CURRENT_SOURCE_DIR}/tDracoInfo.py
  APP    $<TARGET_FILE_DIR:Exe_draco_info>/$<TARGET_FILE_NAME:Exe_draco_info>
  TEST_ARGS "none;--version;--brief"
  LABELS nomemcheck
)

# ---------------------------------------------------------------------------- #
# End diagnostics/test/CMakeLists.txt
# ---------------------------------------------------------------------------- #
