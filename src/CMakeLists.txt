#-----------------------------*-cmake-*----------------------------------------#
# file   src/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 April 28
# brief  Instructions for building src level Makefile.
# note   Copyright (C) 2010-2012 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

cmake_minimum_required( VERSION 2.8 FATAL_ERROR )
project( draco_src_dir )

# Platform Checks 
# Is HOST_NAME_MAX defined?  Is WinSock2.h available?  Is gethostname() available?
message("
Platform Checks...
")
include( platform_checks )
set_draco_uname()
query_have_gethostname()
query_have_maxpathlen()
if( DRACO_ENABLE_CXX11 )
   set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
      "${Draco_SOURCE_DIR}/config/cxx11")
   include(CheckCXX11Features)
   # See ${CXX11_FEATURE_LIST}={HAS_CXX11_AUTO;HAS_CXX11_NULLPTR;HAS_CXX11_LAMBDA;HAS_CXX11_STAT...}
endif()

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
add_dir_if_exists( device )       # needs ds++
add_dir_if_exists( fpe_trap )     # needs ds++
if (NOT APPLE)
  add_dir_if_exists( lapack_wrap )  # needs ds++
endif()
add_dir_if_exists( linear )       # needs ds++
add_dir_if_exists( mesh_element ) # needs ds++
add_dir_if_exists( ode )          # needs ds++
add_dir_if_exists( rng )          # needs ds++
add_dir_if_exists( traits )       # needs ds++
add_dir_if_exists( units )        # needs ds++
if( _LANGUAGES_ MATCHES Fortran )
   if( NOT "${CMAKE_CXX_COMPILER}" MATCHES "[sp]pu-g[+][+]" 
   AND NOT ${CMAKE_GENERATOR} MATCHES Xcode)
      add_dir_if_exists( FortranChecks )
   endif()
endif()

if( GRACE_FOUND )
  add_dir_if_exists( plot2D )
endif()

# Level 3
message(" ")
message( STATUS "Configuring Level 3 packages --" )
add_dir_if_exists( diagnostics )  # needs c4
add_dir_if_exists( fit )          # needs linear
add_dir_if_exists( meshReaders )  # needs c4
add_dir_if_exists( min )          # needs linear
add_dir_if_exists( norms )        # needs c4
add_dir_if_exists( parser )       # needs units, c4
add_dir_if_exists( roots )        # needs linear
add_dir_if_exists( timestep )     # needs c4, ds++
add_dir_if_exists( viz )          # needs traits, ds++
add_dir_if_exists( cdi_ipcress )

# Level 4
message(" ")
message( STATUS "Configuring Level 4 packages --" )
add_dir_if_exists( special_functions ) # needs roots, units, ode
add_dir_if_exists( cdi_analytic )      # needs parser, ode, cdi
#add_dir_if_exists( cdi_eospac )        # needs parser, ode, cdi
add_dir_if_exists( RTT_Format_Reader)  # needs meshReaders

# Level 5
message(" ")
message( STATUS "Configuring Level 5 packages --" )
add_dir_if_exists( quadrature )       # needs mesh_element, parser, special_functions

# Summary
message(" ")

add_feature_info( BUILD_AUTODOC BUILD_AUTODOC
   "Turn on to enable generation of doxygen documentation (make autodoc).")

feature_summary( 
   WHAT PACKAGES_FOUND 
   INCLUDE_QUIET_PACKAGES
   DESCRIPTION "Packages found:")

feature_summary( 
   WHAT PACKAGES_NOT_FOUND 
   INCLUDE_QUIET_PACKAGES
   DESCRIPTION "Packages missing:")

feature_summary( 
   WHAT ENABLED_FEATURES 
   INCLUDE_QUIET_PACKAGES
   DESCRIPTION "Enabled features:")

feature_summary( 
   WHAT DISABLED_FEATURES 
   INCLUDE_QUIET_PACKAGES
   DESCRIPTION "Disabled features:")

message("Draco build summary:

Version     : Draco ${DRACO_VERSION_FULL}
Build type  : ${CMAKE_BUILD_TYPE}
C4 Model    : ${DRACO_C4}
Prefix dir  : ${CMAKE_INSTALL_PREFIX}
Source dir  : ${Draco_SOURCE_DIR}
Build dir   : ${Draco_BINARY_DIR}

CXX Compiler: ${CMAKE_CXX_COMPILER}
CXX FLAGS   : ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}")
if( _LANGUAGES_ MATCHES Fortran) 
   message("FC FLAGS    : ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}}")
endif()
message(
"mpirun cmd  : ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} N ${MPIEXEC_POSTFLAGS}")
message("Library Type: ${DRACO_LIBRARY_TYPE}

ENABLE_RNG_NR    : ${ENABLE_RNG_NR}
DRACO_DBC_LEVEL  : ${DRACO_DBC_LEVEL}")
message("DRACO_DIAGNOSTICS: ${DRACO_DIAGNOSTICS}
DRACO_TIMING     : ${DRACO_TIMING}")

# Push some variables up one level
set( Draco_EXPORT_TARGET_PROPERTIES
   "${Draco_EXPORT_TARGET_PROPERTIES}" PARENT_SCOPE)

##---------------------------------------------------------------------------##
## End of src/CMakeLists.txt
##---------------------------------------------------------------------------##
