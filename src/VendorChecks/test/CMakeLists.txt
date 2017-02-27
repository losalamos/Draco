#-----------------------------*-cmake-*----------------------------------------#
# file   src/VendorsChecks/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2016 May 11
# brief  Instructions for building VendorsChecks/test level Makefile.
# note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
project( VendorsChecks_test CXX )

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories(
   ${PROJECT_SOURCE_DIR}      # headers for tests
   ${Draco_SOURCE_DIR}/src    # headers for package
   ${Draco_BINARY_DIR}/src/c4      # c4/config.h
   ${Draco_BINARY_DIR}/src/ds++ )  # ds++/config.h

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

if( METIS_FOUND )
  add_scalar_tests(
    SOURCES "tstMetis.cc"
    DEPS    "Lib_dsxx;METIS::metis"
    )
endif()

# Only run these tests if MPI is found and enabled. Both parmetis and
# superlu-dist require MPI.
if( ${DRACO_C4} STREQUAL "MPI" )

  if( ParMETIS_FOUND )
    add_parallel_tests(
      SOURCES "tstParmetis.cc"
      DEPS    "Lib_c4;ParMETIS::parmetis"
      PE_LIST "3"
      )
  endif()

  if( SuperLU_DIST_FOUND )
    add_parallel_tests(
      SOURCES "tstSuperludist.cc"
      DEPS    "Lib_c4;SuperLU_DIST::superludist;ParMETIS::parmetis;lapack"
      PE_LIST "4"
      )
  endif()

endif()

# ---------------------------------------------------------------------------- #
# end src/VendorsChecks/test/CMakeLists.txt
# ---------------------------------------------------------------------------- #
