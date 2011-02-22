#=============================================================================
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
#=============================================================================
# $Id$
#=============================================================================

project( parser_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB test_lib_sources *_pt.cc )
list( APPEND test_lib_sources 
   ${PROJECT_SOURCE_DIR}/parser_test.cc 
   ${PROJECT_SOURCE_DIR}/parser_test.hh)
file( GLOB test_sources *.cc )
list( REMOVE_ITEM test_sources ${test_lib_sources} )
file( GLOB input_files *.inp )

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories( 
   ${PROJECT_SOURCE_DIR}      # headers for tests
   ${PROJECT_SOURCE_DIR}/..   # headers for package
   ${PROJECT_BINARY_DIR}/..   # config.h
   ${MPI_INCLUDE_PATH} )

# Library Lib_parser_test
add_component_library( Lib_parser_test ${PROJECT_NAME} "${test_lib_sources}" )
add_dependencies( Lib_parser_test Lib_parser )

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
     PROJECT_LABEL tstConsole_Token_Stream
   )
target_link_libraries( Ut_parser_tstConsole_Token_Stream
   Lib_parser_test 
   Lib_parser
   Lib_c4
   Lib_dsxx
   ${MPI_LIBRARIES}
   ${INTEL_LIBM}
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
   DEPS    "Lib_c4;Lib_dsxx;${MPI_LIBRARIES};${INTEL_LIBM}" )

# Special case for driver4tstConsole_Token_Stream
# 1. This test doesn't work under Xcode
# 2. This test must have unique access to tstConsole_Token_Stream and
#    we will use a RESOURCE_LOCK to prevent multiple simultaneious
#    instances during parallel ctest.
list( REMOVE_ITEM test_sources 
   ${PROJECT_SOURCE_DIR}/driver4tstConsole_Token_Stream.cc )
if( NOT ${CMAKE_GENERATOR} MATCHES "Xcode" )
add_parallel_tests(
   SOURCES "${PROJECT_SOURCE_DIR}/driver4tstConsole_Token_Stream.cc"
   PE_LIST "1;2"
   DEPS    "Lib_c4;Lib_dsxx;${MPI_LIBRARIES}" 
   RESOURCE_LOCK "lock_tstConsole_Token_Stream" )
endif()

# The remaining tests are Parallel.

add_parallel_tests(
   SOURCES "${test_sources}"
   PE_LIST "1;2"
   DEPS    "Lib_c4;Lib_dsxx;${MPI_LIBRARIES};${INTEL_LIBM}" )
   
# ---------------------------------------------------------------------------- #
# Copy input files to test directory (binary_dir)
# ---------------------------------------------------------------------------- #

provide_aux_files( FILES "${input_files}" )
   
