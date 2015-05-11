#-----------------------------*-cmake-*----------------------------------------#
# file   parser/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for parser/test.
# note   Copyright (C) 2010-2015, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#
project( parser_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB test_sources *.cc )
file( GLOB input_files *.inp )

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

# Special Case:
# tstConsole_Token_Stream is used by driver4tstConsole_Token_Stream,
# we do not run it alone.

add_executable( Ut_parser_tstConsole_Token_Stream
   ${PROJECT_SOURCE_DIR}/tstConsole_Token_Stream.cc )
set_target_properties( Ut_parser_tstConsole_Token_Stream
   PROPERTIES
     OUTPUT_NAME tstConsole_Token_Stream
     VS_KEYWORD  tstConsole_Token_Stream
     PROJECT_LABEL Ut_parser_tstConsole_Token_Stream
     FOLDER      parser_test
   )
target_link_libraries( Ut_parser_tstConsole_Token_Stream
   Lib_parser
   Lib_c4
   Lib_dsxx
   ${MPI_CXX_LIBRARIES}
   ${PAPI_LIBRARY}
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
list( REMOVE_ITEM test_sources
   ${PROJECT_SOURCE_DIR}/driver4tstConsole_Token_Stream.cc )
if( NOT ${CMAKE_GENERATOR} MATCHES "Xcode" )
add_scalar_tests(
   SOURCES "${PROJECT_SOURCE_DIR}/driver4tstConsole_Token_Stream.cc"
   DEPS    "Lib_parser;${MPI_CXX_LIBRARIES};${PAPI_LIBRARY}"
   RESOURCE_LOCK "lock_tstConsole_Token_Stream" )
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

# ---------------------------------------------------------------------------- #
# Copy input files to test directory (binary_dir)
# ---------------------------------------------------------------------------- #

provide_aux_files( FILES "${input_files}" )
