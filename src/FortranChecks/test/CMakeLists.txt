# ---------------------------------------------------------------------------- #
# file   FortranChecks/test/CMakeLists.txt
# author Kelly Thompson
# brief  Instructions for building a C++ main that links against an
#        F90 library and uses ISO_C_BINDING.
# note   Copyright (C) 2016-2019 Triad National Security, LLC.
#        All rights reserved
#------------------------------------------------------------------------------#
project( FortranChecks_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #
set( test_sources
  ${PROJECT_SOURCE_DIR}/cppmain.cc
  ${PROJECT_SOURCE_DIR}/tstDerivedTypes.cc
  ${PROJECT_SOURCE_DIR}/tstDrel.cc
  )
set( ptest_sources
  ${PROJECT_SOURCE_DIR}/mpi_hw.cc
)

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #
# include_directories(
#   ${PROJECT_SOURCE_DIR}      # headers for tests
#   ${Draco_SOURCE_DIR}/src    # headers for ds++
#   )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

# When building on SQ with XLC/XLF90, ensure that libxlomp_ser is not
# on the link line.  This macro is defined in config/unix-xlf.cmake.
if( ${CMAKE_Fortran_COMPILER_FLAVOR} MATCHES "XL" )
  remove_lib_from_link( "xlomp_ser" )
endif()

# When building with MVSE and linking Fortran libraries generated with MinGW's
# gfortran, the libquadmath-0.dll library must be found in the path.
if( MSVC )
  find_file( mingw_libquadmath libquadmath-0.dll )
  if( NOT mingw_libquadmath )
    message( FATAL_ERROR "ERROR in FortranChecks\test:
      The library libquadmath-0.dll must be found in your PATH to allow these tests to run!")
  endif()
endif()

# Create the test.
set( test_deps "Lib_FC_f90sub" )

#include(print_target_properties)
#print_targets_properties( "Lib_FC_f90sub;Lib_FC_Derived_Type;Lib_dsxx" )

# When cross-compiling with '-mmic', rpaths for libraries built from
# Fortran code don't appear to be reported to the icpc linker.  As a
# hack, these libraries (discovered in unix-ifort.cmake) are manually
# added to the link line.
if( HAVE_MIC )
  list( APPEND test_deps ${libifport_loc} ${libifcoremt_loc} )
endif()

add_parallel_tests(
  SOURCES "${test_sources}"
  DEPS    "${test_deps}"
  PE_LIST "1" )

add_parallel_tests(
  SOURCES "${ptest_sources}"
  DEPS    "${test_deps}"
  PE_LIST "1;2" )

# ---------------------------------------------------------------------------- #
# end of FortranChecks/test/CMakeLists.txt
# ---------------------------------------------------------------------------- #
