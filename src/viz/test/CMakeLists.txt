#-----------------------------*-cmake-*----------------------------------------#
# file   viz/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for viz/test.
# note   Copyright (C) 2010-2015, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#
project( viz_test CXX )

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

if( NOT APPLE )
   # numdiff/ndselect crashes when run on OS/X so don't try to run the
   # tests.
   add_scalar_tests(
      SOURCES "${test_sources}"
      DEPS    "Lib_viz;Lib_dsxx"
      )
endif()

##---------------------------------------------------------------------------##
## Extra 'make clean' files
##---------------------------------------------------------------------------##
set( extra_clean_files
   Densities.diff
   Pressure.diff
   Temperatures.diff
   Velocity.diff
   ensight_stream.out
   geo.diff
   part_testproblem_binary_ensight
   part_testproblem_ensight
   testproblem_binary_ensight
   testproblem_ensight
)
set_directory_properties(
   PROPERTIES
   ADDITIONAL_MAKE_CLEAN_FILES "${extra_clean_files}" )
