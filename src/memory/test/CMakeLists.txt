#-----------------------------*-cmake-*----------------------------------------#
# file   config/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for memory/test.
# note   Copyright (C) 2010-2014, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id: CMakeLists.txt 7204 2013-09-03 14:38:22Z kellyt $
#------------------------------------------------------------------------------#
project( memory_test CXX )

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
   )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_scalar_tests(
   SOURCES "${test_sources}" 
   DEPS    Lib_memory
   )
   
# Application Unit Test 
add_test( 
   NAME    memory_draco_info
   COMMAND ${CMAKE_COMMAND} 
   -DAPP=$<TARGET_FILE_DIR:Exe_draco_info>/$<TARGET_FILE_NAME:Exe_draco_info>
   -DWORKDIR=$<TARGET_FILE_DIR:Ut_memory_tstmemory_exe>
   -P ${PROJECT_SOURCE_DIR}/tDracoInfo.cmake
   )
set_tests_properties( memory_draco_info
   PROPERTIES	
     PASS_REGULAR_EXPRESSION "Passes"
     FAIL_REGULAR_EXPRESSION "Fails"
     LABELS "nomemcheck" # valgrind gets confused when
                         # parsing a 32-bit program (e.g.: cmake)
   )

