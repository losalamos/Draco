#-----------------------------*-cmake-*----------------------------------------#
# file   src/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 April 28
# brief  Instructions for building src level Makefile.
# note   Copyright (C) 2016-2019 Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.9.0)
project( draco_src_dir )

# Provide helper functions used by component CMakeLists.txt files
include( component_macros )

# Extra 'draco-only' flags (don't pass these to Capsaicin).
if( CMAKE_CXX_COMPILER_ID STREQUAL "GNU" )
  if( CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 5.0 )
    toggle_compiler_flag( TRUE "-Wfloat-equal -Wunused-macros" "CXX;C" "DEBUG")
    # toggle_compiler_flag( TRUE "-Winline" "CXX;C" "RELEASE")
  endif()
elseif( CMAKE_CXX_COMPILER_ID STREQUAL "MSVC" )
  string( REPLACE "/W2" "/W4" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}" )
  string( REPLACE "/W2" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" )
endif()

# Level 1
message(" ")
message( STATUS "Configuring Level 1 packages --" )
add_dir_if_exists( ds++ )

# Level 2
message(" ")
message( STATUS "Configuring Level 2 packages --" )
add_dir_if_exists( c4 )           # needs ds++
add_dir_if_exists( cdi )          # needs ds++
add_dir_if_exists( linear )       # needs ds++
add_dir_if_exists( mesh_element ) # needs ds++
add_dir_if_exists( ode )          # needs ds++
add_dir_if_exists( rng )          # needs ds++
add_dir_if_exists( units )        # needs ds++
add_dir_if_exists( viz )          # needs ds++
add_dir_if_exists( memory )       # needs ds++

# For now, do not try to build or test some packages (these two paths are
# identical now, but are expected to diverge as more packages are added to the
# CI test system):
#  - For appveyor (Visual Studio) CI on Github, omit packages that need special
#    libraries like CSK, lapack, metis, parmetis or cuda. Also omit packages
#    that require Fortran.
if( NOT DEFINED ENV{APPVEYOR} )
  add_dir_if_exists( device )       # needs ds++
  add_dir_if_exists( compton )      # needs c4, CSK
  if( HAVE_Fortran )
    add_dir_if_exists( lapack_wrap )  # needs ds++
  endif()
  if( Grace_FOUND )
    add_dir_if_exists( plot2D )
  endif()
endif()

# Level 3
message(" ")
message( STATUS "Configuring Level 3 packages --" )
add_dir_if_exists( cdi_ipcress )
add_dir_if_exists( diagnostics )  # needs c4
add_dir_if_exists( fit )          # needs linear
add_dir_if_exists( meshReaders )  # needs c4
add_dir_if_exists( min )          # needs linear
add_dir_if_exists( norms )        # needs c4
add_dir_if_exists( parser )       # needs units, c4
add_dir_if_exists( roots )        # needs linear
add_dir_if_exists( timestep )     # needs c4, ds++

if( NOT DEFINED ENV{APPVEYOR} )
  add_dir_if_exists( VendorChecks ) # needs c4
  if( HAVE_Fortran )
    add_dir_if_exists( FortranChecks )
  endif()
  if( NOT DEFINED ENV{TRAVIS} )
    add_dir_if_exists( cdi_eospac ) # needs parser, ode, cdi
  endif()
endif()

# Level 4
message(" ")
message( STATUS "Configuring Level 4 packages --" )
add_dir_if_exists( special_functions ) # needs roots, units, GSL
add_dir_if_exists( cdi_analytic )      # needs parser, roots, ode, cdi
add_dir_if_exists( RTT_Format_Reader)  # needs meshReaders

# Level 5
message(" ")
message( STATUS "Configuring Level 5 packages --" )
if( NOT DEFINED ENV{APPVEYOR} )
  add_dir_if_exists( quadrature ) # needs mesh_element, parser, special_functions
endif()
add_dir_if_exists( mesh )         # needs c4, ds++, mesh_element

# Summary

message(" ")
feature_summary( WHAT ALL
  INCLUDE_QUIET_PACKAGES
  FATAL_ON_MISSING_REQUIRED_PACKAGES
  QUIET_ON_EMPTY )

message("
Draco build summary:

Version     : Draco ${Draco_VERSION_FULL}
Build type  : ${Draco_BUILD_TYPE}
            : DRACO_DBC_LEVEL   = ${DRACO_DBC_LEVEL}
            : DRACO_DIAGNOSTICS = ${DRACO_DIAGNOSTICS}
            : DRACO_TIMING      = ${DRACO_TIMING}
C4 Model    : ${DRACO_C4}
Prefix dir  : ${CMAKE_INSTALL_PREFIX}
Source dir  : ${Draco_SOURCE_DIR}
Build dir   : ${Draco_BINARY_DIR}

CXX Compiler: ${CMAKE_CXX_COMPILER}")
if( CMAKE_CONFIGURATION_TYPES )
  message("CXX FLAGS   : ${CMAKE_CXX_FLAGS} ")
  message("CXX Debug FL: ${CMAKE_CXX_FLAGS_DEBUG}")
  message("CXX Release : ${CMAKE_CXX_FLAGS_RELEASE}")
  if( _LANGUAGES_ MATCHES Fortran)
    message("Fortran     : ${CMAKE_Fortran_COMPILER}")
    message("FC FLAGS    : ${CMAKE_Fortran_FLAGS}")
    message("FC Debug FLA: ${CMAKE_Fortran_FLAGS_DEBUG}")
    message("FC Release F: ${CMAKE_Fortran_FLAGS_RELEASE}")
  endif()
else()
  message("C FLAGS     : ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${Draco_BUILD_TYPE}}")
  message("CXX FLAGS   : ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${Draco_BUILD_TYPE}}")
  if( _LANGUAGES_ MATCHES Fortran)
    message("Fortran     : ${CMAKE_Fortran_COMPILER}")
    message("FC FLAGS    : ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${Draco_BUILD_TYPE}}")
  endif()
endif()
if( CAFS_Fortran_COMPILER )
  message("CAFS Fortran: ${CAFS_Fortran_COMPILER}")
endif()
if( "${DRACO_C4}" STREQUAL "MPI" )
  message(
  "mpirun cmd  : ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} N ${MPIEXEC_POSTFLAGS}")
endif()
message("Library Type: ${DRACO_LIBRARY_TYPE}
")
if( CRAY_PE AND ENV{CRAYPE_VERSION} )
  message("Cray system detected: CC -craype-verbose -V:
")
  if( ${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel" )
    set( ver_opt "--version")
  else( ${CMAKE_CXX_COMPILER_ID} STREQUAL "Cray" )
    set( ver_opt "-V")
  endif()
  execute_process(
    COMMAND ${CMAKE_CXX_COMPILER} -craype-verbose ${ver_opt}
    OUTPUT_VARIABLE tmp
    )
  message("${tmp}")
endif()

# Push some variables up one level
set( Draco_EXPORT_TARGET_PROPERTIES
  "${Draco_EXPORT_TARGET_PROPERTIES}" PARENT_SCOPE)

##---------------------------------------------------------------------------##
## Debug build system state:
##---------------------------------------------------------------------------##
option( DBS_PRINT_DEBUG_INFO "Print extra build system data" OFF )
if( DBS_PRINT_DEBUG_INFO )
  include( debug_macros )
  # Debug variables
  message( "CMake variables and values:" )
  echo_all_cmake_variable_values()

  # Debug target properties (run 'make help' to see a list of targets)
  message( "CMake targets and properties:" )
  set( dbs_targets
    # Lib_dsxx
    # Lib_dsxx_test
    Exe_Ipcress_Interpreter
    )
  echo_targets( ${dbs_targets} )
endif()

##---------------------------------------------------------------------------##
## End of src/CMakeLists.txt
##---------------------------------------------------------------------------##
