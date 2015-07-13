#-----------------------------*-cmake-*----------------------------------------#
# file   cdi_eospac/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for cdi_eospac/test.
# note   Copyright (C) 2010-2015, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id: CMakeLists.txt 6732 2012-09-05 22:28:18Z kellyt $
#------------------------------------------------------------------------------#
project( cdi_eospac_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_sources
   ${PROJECT_SOURCE_DIR}/tEospac.cc
   ${PROJECT_SOURCE_DIR}/tEospacWithCDI.cc
)

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories(
   ${PROJECT_SOURCE_DIR}      # headers for tests
   ${PROJECT_SOURCE_DIR}/..   # headers for package
   ${PROJECT_BINARY_DIR}/..   # c4/config.h
)

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_scalar_tests(
   SOURCES "${test_sources}"
   DEPS    "Lib_cdi_eospac"
   )

add_test(
   NAME cdi_eospac_tQueryEospac
   COMMAND ${CMAKE_COMMAND}
   -DAPP=$<TARGET_FILE_DIR:Exe_QueryEospac>/$<TARGET_FILE_NAME:Exe_QueryEospac>
   -DSTDINFILE=${PROJECT_SOURCE_DIR}/QueryEospacInput.dat
   -DGOLDFILE=${PROJECT_SOURCE_DIR}/QueryEospac.gold
   -DDraco_BINARY_DIR=${Draco_BINARY_DIR}
   -P ${PROJECT_SOURCE_DIR}/tQueryEospac.cmake
)
set_tests_properties( cdi_eospac_tQueryEospac
   PROPERTIES
      PASS_REGULAR_EXPRESSION "Passes"
      FAIL_REGULAR_EXPRESSION "Fails"
      LABELS "nomemcheck"  # valgrind gets confused when
                           # parsing a 32-bit program (e.g.: cmake)
)

# Test QueryEospac '--version' option
add_test(
   NAME cdi_eospac_tQueryEospac_version
   COMMAND ${CMAKE_COMMAND}
   -DAPP=$<TARGET_FILE_DIR:Exe_QueryEospac>/$<TARGET_FILE_NAME:Exe_QueryEospac>
   -DCMDARGS=--version
   -DDraco_BINARY_DIR=${Draco_BINARY_DIR}
   -P ${PROJECT_SOURCE_DIR}/tQueryEospac.cmake
)
set_tests_properties( cdi_eospac_tQueryEospac_version
   PROPERTIES
      PASS_REGULAR_EXPRESSION "Passes"
      FAIL_REGULAR_EXPRESSION "Fails"
)

# Test QueryEospac '--help' option
add_test(
   NAME cdi_eospac_tQueryEospac_help
   COMMAND ${CMAKE_COMMAND}
   -DAPP=$<TARGET_FILE_DIR:Exe_QueryEospac>/$<TARGET_FILE_NAME:Exe_QueryEospac>
   -DCMDARGS=--help
   -DDraco_BINARY_DIR=${Draco_BINARY_DIR}
   -P ${PROJECT_SOURCE_DIR}/tQueryEospac.cmake
)
set_tests_properties( cdi_eospac_tQueryEospac_help
   PROPERTIES
      PASS_REGULAR_EXPRESSION "Passes"
      FAIL_REGULAR_EXPRESSION "Fails"
)
