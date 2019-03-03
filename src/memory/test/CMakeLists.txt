#-----------------------------*-cmake-*----------------------------------------#
# file   memory/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for memory/test.
# note   Copyright (C) 2016-2019, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
project( memory_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB test_sources *.cc )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_scalar_tests( SOURCES "${test_sources}"  DEPS Lib_memory  )

#------------------------------------------------------------------------------#
# End memory/test/CMakeLists.txt
#------------------------------------------------------------------------------#
