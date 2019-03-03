#-----------------------------*-cmake-*----------------------------------------#
# file   norms/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for norms/test.
# note   Copyright (C) 2016-2019, Triad National Security
#        All rights reserved.
#------------------------------------------------------------------------------#
project( norms_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB test_sources *.cc )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_parallel_tests(
   SOURCES "${test_sources}"
   PE_LIST "1;2"
#   DEPS    "Lib_norms;${MPI_CXX_LIBRARIES};${PAPI_LIBRARY}" )
   DEPS    "Lib_norms;${PAPI_LIBRARY}" )

#------------------------------------------------------------------------------#
# End norms/test/CMakeLists.txt
#------------------------------------------------------------------------------#
