#-----------------------------*-cmake-*----------------------------------------#
# file   src/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 April 28
# brief  Instructions for building src level Makefile.
# note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.0.0)
project( draco_src_dir )

# Provide helper functions used by component CMakeLists.txt files
include( component_macros )

# Helper macro that only tries to configure a listed package if it
# exists in the source tree.  For example, when building milagro, the
# get_draco script only checks out the parts of Draco that are
# required by milagro.
macro( add_dir_if_exists package )
  message( "   ${package}" )
  if( EXISTS ${PROJECT_SOURCE_DIR}/${package} )
    add_subdirectory( ${package} )
  endif()
endmacro()

# Level 1
message(" ")
message( STATUS "Configuring Level 1 packages --" )
add_dir_if_exists( ds++ )

# Level 2
message(" ")
message( STATUS "Configuring Level 2 packages --" )
add_dir_if_exists( c4 )           # needs ds++
add_dir_if_exists( cdi )          # needs ds++
add_dir_if_exists( compton )      # needs ds++
add_dir_if_exists( device )       # needs ds++
add_dir_if_exists( linear )       # needs ds++
add_dir_if_exists( mesh_element ) # needs ds++
add_dir_if_exists( ode )          # needs ds++
add_dir_if_exists( rng )          # needs ds++
add_dir_if_exists( units )        # needs ds++
add_dir_if_exists( viz )          # needs ds++
if( HAVE_Fortran AND
    ${CMAKE_GENERATOR} MATCHES "Unix Makefiles" AND
    NOT CMAKE_C_COMPILER MATCHES clang )
  add_dir_if_exists( lapack_wrap )  # needs ds++
endif()
if( GRACE_FOUND )
  add_dir_if_exists( plot2D )
endif()
if( NOT WIN32 )
  add_dir_if_exists( memory )         # needs dsxx
endif()

# Level 3
message(" ")
message( STATUS "Configuring Level 3 packages --" )
add_dir_if_exists( VendorChecks ) # needs c4
add_dir_if_exists( cdi_ipcress )
add_dir_if_exists( diagnostics )  # needs c4
add_dir_if_exists( fit )          # needs linear
add_dir_if_exists( meshReaders )  # needs c4
add_dir_if_exists( min )          # needs linear
add_dir_if_exists( norms )        # needs c4
add_dir_if_exists( parser )       # needs units, c4
add_dir_if_exists( roots )        # needs linear
add_dir_if_exists( timestep )     # needs c4, ds++
if( HAVE_Fortran )
  add_dir_if_exists( FortranChecks )
endif()

# Level 4
message(" ")
message( STATUS "Configuring Level 4 packages --" )
add_dir_if_exists( special_functions ) # needs roots, units, ode
add_dir_if_exists( cdi_analytic )      # needs parser, ode, cdi
add_dir_if_exists( cdi_eospac )        # needs parser, ode, cdi
add_dir_if_exists( RTT_Format_Reader)  # needs meshReaders

# Level 5
message(" ")
message( STATUS "Configuring Level 5 packages --" )
add_dir_if_exists( quadrature )       # needs mesh_element, parser, special_functions

# Summary
message(" ")

feature_summary(
  WHAT PACKAGES_FOUND
  INCLUDE_QUIET_PACKAGES
  DESCRIPTION "Packages found:")

feature_summary(
  WHAT PACKAGES_NOT_FOUND
  INCLUDE_QUIET_PACKAGES
  FATAL_ON_MISSING_REQUIRED_PACKAGES
  DESCRIPTION "Packages missing:")

feature_summary(
  WHAT ENABLED_FEATURES
  INCLUDE_QUIET_PACKAGES
  DESCRIPTION "Enabled features:")

feature_summary(
  WHAT DISABLED_FEATURES
  INCLUDE_QUIET_PACKAGES
#  QUIET_ON_EMPTY # only available for cmake-3.8+
  DESCRIPTION "Disabled features:")

message("
Draco build summary:

Version     : Draco ${DRACO_VERSION_FULL}
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
message(
  "mpirun cmd  : ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} N ${MPIEXEC_POSTFLAGS}")
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
