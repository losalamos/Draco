#-----------------------------*-cmake-*----------------------------------------#
# file   src/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 April 28
# brief  Instructions for building src level Makefile.
# note   © Copyright 2010 Los Alamos National Security, All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

cmake_minimum_required( VERSION 2.8 FATAL_ERROR )
project( draco_src_dir )

# Set default compile environment:

# Setup defaults, value checks, etc.
include(buildEnv)
# set defaults for BUILD_TYPE and INSTALL_PREFIX
dbsSetDefaults()
# Save config info
dbsConfigInfo()

# Set compiler options
include( compilerEnv )
dbsSetupCompilers()
dbsSetupCxx()

# Find any globally required libraries
include( vendor_libraries )
setupVendorLibraries()

#
# Extra Make targets (package_commands.cmake)
#
#provide_test_target() 

# Level 1
message( STATUS "\nConfiguring Level 1 packages --\n" )
add_subdirectory( ds++ )

# Level 2
message( STATUS "\nConfiguring Level 2 packages --\n" )
add_subdirectory( c4 )           # needs ds++
add_subdirectory( cdi )          # needs ds++
add_subdirectory( fpe_trap )     # needs ds++
add_subdirectory( lapack_wrap )  # needs ds++
add_subdirectory( linear )       # needs ds++
add_subdirectory( mesh_element ) # needs ds++
add_subdirectory( ode )          # needs ds++
add_subdirectory( rng )          # needs ds++
add_subdirectory( traits )       # needs ds++
add_subdirectory( units )        # needs ds++
add_subdirectory( xm )           # needs ds++

if( PCG_FOUND AND EXISTS ${PROJECT_SOURCE_DIR}/pcgWrap/CMakeLists.txt)
  add_subdirectory( pcgWrap )    # needs ds++
endif()
if( XMGRACE_FOUND AND EXISTS ${PROJECT_SOURCE_DIR}/plot2D/CMakeLists.txt)
  add_subdirectory( plot2D )
endif()

# Level 3

message( STATUS "\nConfiguring Level 3 packages --\n" )
add_subdirectory( cdi_analytic ) # needs cdi
add_subdirectory( cdi_gandolf )  # needs cdi
add_subdirectory( diagnostics )  # needs c4
add_subdirectory( fit )          # needs lienar
add_subdirectory( meshReaders )  # needs c4
add_subdirectory( min )          # needs linear
add_subdirectory( norms )        # needs c4
add_subdirectory( parser )       # needs units, c4
add_subdirectory( roots  )       # needs linear
add_subdirectory( special_functions ) # needs units, ds++
add_subdirectory( timestep )     # needs c4, ds++
add_subdirectory( viz )          # needs traits, ds++

# Level 4

message( STATUS "\nConfiguring Level 4 packages --\n" )
add_subdirectory( quadrature ) # needs units, c4, special_functions
add_subdirectory( RTT_Format_Reader) # needs meshReaders, mesh_element, c4, dsxx

# Summary
message("
Draco build summary:

Version     : Draco ${DRACO_VERSION_FULL}
Build type  : ${CMAKE_BUILD_TYPE}
C4 Model    : ${DRACO_C4}
Prefix dir  : ${CMAKE_INSTALL_PREFIX}
Source dir  : ${Draco_SOURCE_DIR}
Build dir   : ${Draco_BINARY_DIR}

CXX Compiler: ${CMAKE_CXX_COMPILER}
CXX FLAGS   : ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}
Library Type: ${DRACO_LIBRARY_TYPE}

Vendor dir  : ${VENDOR_DIR}
MPI         : ${MPIEXEC}
GSL         : ${GSL_LIBRARY}
LAPACK      : ${LAPACK_LIBRARY}
GANDOLF     : ${GANDOLF_LIBRARY}

")

## Install cmake files that help external projects link to the targets from this project.
#install( 
#  EXPORT draco_targets
#  NAMESPACE draco_
#  DESTINATION lib/export
#)

## Build summary
#if( VERBOSE )
#  include(summary)
#endif( VERBOSE )

##---------------------------------------------------------------------------##
## End of src/CMakeLists.txt
##---------------------------------------------------------------------------##
