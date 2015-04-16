#-----------------------------*-cmake-*----------------------------------------#
# file   config/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for diagnostics/test.
# note   Copyright (C) 2010-2015, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
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

# Application Unit Test
add_test(
   NAME    diagnostics_draco_info
   COMMAND ${CMAKE_COMMAND}
   -DAPP=$<TARGET_FILE_DIR:Exe_draco_info>/$<TARGET_FILE_NAME:Exe_draco_info>
   -DWORKDIR=$<TARGET_FILE_DIR:Ut_diagnostics_tstDiagnostics_exe>
   -DHAVE_MIC=${HAVE_MIC}
   -DDraco_BINARY_DIR=${Draco_BINARY_DIR}
   -P ${Draco_SOURCE_DIR}/src/diagnostics/test/tDracoInfo.cmake
   )
set_tests_properties( diagnostics_draco_info
   PROPERTIES
     PASS_REGULAR_EXPRESSION "Passes"
     FAIL_REGULAR_EXPRESSION "Fails"
     LABELS "nomemcheck" # valgrind gets confused when
                         # parsing a 32-bit program (e.g.: cmake)
   )
