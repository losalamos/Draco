#-----------------------------*-cmake-*----------------------------------------#
# file   parser/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for parser/test.
# note   Copyright (C) 2016-2017, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
project( parser_test CXX )

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

add_component_executable(
  TARGET      Exe_tstConsole_Token_Stream
  SOURCES     ${PROJECT_SOURCE_DIR}/tstConsole_Token_Stream.cc
  TARGET_DEPS Lib_parser
  NOEXPORT
  )
list( REMOVE_ITEM test_sources
   ${PROJECT_SOURCE_DIR}/tstConsole_Token_Stream.cc )

# These tests are scalar only
set( scalar_test_sources
   ${PROJECT_SOURCE_DIR}/tstToken.cc
   ${PROJECT_SOURCE_DIR}/tstToken_Equivalence.cc
   ${PROJECT_SOURCE_DIR}/tstutilities.cc
   ${PROJECT_SOURCE_DIR}/tstUnit.cc
   ${PROJECT_SOURCE_DIR}/tstString_Token_Stream.cc
   ${PROJECT_SOURCE_DIR}/tstParse_Table.cc )
list( REMOVE_ITEM test_sources ${scalar_test_sources} )
add_scalar_tests(
   SOURCES "${scalar_test_sources}"
   DEPS    "Lib_parser;${MPI_CXX_LIBRARIES};${PAPI_LIBRARY}" )

# Special case for driver4tstConsole_Token_Stream
# 1. This test doesn't work under Xcode
# 2. This test must have unique access to tstConsole_Token_Stream and
#    we will use a RESOURCE_LOCK to prevent multiple simultaneous
#    instances during parallel ctest.
if( NOT ${CMAKE_GENERATOR} MATCHES "Xcode" )
  include( ApplicationUnitTest )
  add_app_unit_test(
    DRIVER    ${CMAKE_CURRENT_SOURCE_DIR}/driver4tstConsole_Token_Stream.py
    APP       $<TARGET_FILE_DIR:Exe_tstConsole_Token_Stream>/$<TARGET_FILE_NAME:Exe_tstConsole_Token_Stream>
    STDINFILE "${CMAKE_CURRENT_SOURCE_DIR}/console_test.inp" )
endif()

# The remaining tests are Parallel.

# Run with with 1 PE only
list( REMOVE_ITEM test_sources
   ${PROJECT_SOURCE_DIR}/tstFile_Token_Stream.cc )
add_parallel_tests(
   SOURCES "${PROJECT_SOURCE_DIR}/tstFile_Token_Stream.cc"
   PE_LIST "1"
   DEPS    "Lib_parser;${MPI_CXX_LIBRARIES};${PAPI_LIBRARY}" )

# Run remaining tests with 1 and 2 PEs.
add_parallel_tests(
   SOURCES "${test_sources}"
   PE_LIST "1;2"
   DEPS    "Lib_parser;${MPI_CXX_LIBRARIES};${PAPI_LIBRARY}" )

#------------------------------------------------------------------------------#
# End parser/CMakeLists.txt
#------------------------------------------------------------------------------#
