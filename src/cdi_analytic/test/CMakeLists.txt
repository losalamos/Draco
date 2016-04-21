#-----------------------------*-cmake-*----------------------------------------#
# file   config/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for cdi_analytic/test.
# note   Copyright (C) 2016, Los Alamos National Security
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#
project( cdi_analytic_test CXX )

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
   ${PROJECT_BINARY_DIR}/..   # c4/config.h
)

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #
set ( parallel_test
      ${PROJECT_SOURCE_DIR}/tstAnalytic_Odfmg_Opacity.cc
      ${PROJECT_SOURCE_DIR}/tstPseudo_Line_Analytic_MultigroupOpacity.cc )
add_parallel_tests(
   SOURCES "${parallel_test}"
   DEPS    "Lib_cdi_analytic"
   PE_LIST "1"
)

list ( REMOVE_ITEM test_sources
    ${PROJECT_SOURCE_DIR}/tstAnalytic_Odfmg_Opacity.cc
    ${PROJECT_SOURCE_DIR}/tstPseudo_Line_Analytic_MultigroupOpacity.cc)


add_scalar_tests(
   SOURCES "${test_sources}"
   DEPS    "Lib_cdi_analytic" )
