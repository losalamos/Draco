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

# Build only for Opteron? Or PPE also?
if( "${SITE}" MATCHES "^rr.*" ) # rr-dev-fe, rro*
   include( platform_customization_roadrunner_ppe )
endif()

# Find any globally required libraries
include( FeatureSummary )
include( vendor_libraries )
setupVendorLibraries()

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
add_dir_if_exists( fpe_trap )     # needs ds++
add_dir_if_exists( lapack_wrap )  # needs ds++
add_dir_if_exists( linear )       # needs ds++
add_dir_if_exists( mesh_element ) # needs ds++
add_dir_if_exists( ode )          # needs ds++
add_dir_if_exists( rng )          # needs ds++
add_dir_if_exists( traits )       # needs ds++
add_dir_if_exists( units )        # needs ds++

if( ENABLE_HET_CODE )
   add_dir_if_exists( device )    # needs ds++
endif()
if( PCG_FOUND )
  add_dir_if_exists( pcgWrap )    # needs ds++
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
add_dir_if_exists( special_functions ) # needs units, ds++
add_dir_if_exists( timestep )     # needs c4, ds++
add_dir_if_exists( viz )          # needs traits, ds++

if( GANDOLF_FOUND )
   add_dir_if_exists( cdi_gandolf )  # needs cdi
endif()

# Level 4
message(" ")
message( STATUS "Configuring Level 4 packages --" )
add_dir_if_exists( cdi_analytic )     # needs cdi, ode, parser
add_dir_if_exists( quadrature )       # needs units, c4, special_functions
add_dir_if_exists( RTT_Format_Reader) # needs meshReaders, mesh_element, c4, 
                                      #       dsxx

# Summary
message(" ")

add_feature_info( ENABLE_HET_CODE ENABLE_HET_CODE 
   "Build software for heterogeneous architectures (e.g. roadrunner).
   Required for building the device component.")
add_feature_info( BUILD_AUTODOC BUILD_AUTODOC
   "Turn on to enable generation of doxygen documentation (make autodoc).")
add_feature_info( BUILD_TESTING BUILD_TESTING
   "Turn off to prevent the compilation of unit tests (ctest).")

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
CXX FLAGS   : ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}
Library Type: ${DRACO_LIBRARY_TYPE}

ENABLE_RNG_NR    : ${ENABLE_RNG_NR}
DRACO_DBC_LEVEL  : ${DRACO_DBC_LEVEL}")
#message("ENABLE_HET_CODE  : ${ENABLE_HET_CODE}")
message("DRACO_DIAGNOSTICS: ${DRACO_DIAGNOSTICS}
DRACO_TIMING     : ${DRACO_TIMING}")
# messge("Vendor dir  : ${VENDOR_DIR}
# MPI         : ${MPIEXEC}
# GSL         : ${GSL_LIBRARY}")

# if( LAPACK_FOUND )
#    set(tmplist ${LAPACK_LIBRARIES})
#    list(GET tmplist 1 lib1)
#    message("LAPACK      : ${lib1}")
#    list(LENGTH tmplist list_length)
#    if( ${list_length} GREATER 1 )
#       list(REMOVE_AT tmplist 1 ) 
#       foreach(lib ${tmplist} )
#          message("              ${lib}")
#       endforeach()
#    endif()
# else()
#    message("LAPACK      : NOTFOUND")
# endif()
# message("GANDOLF     : ${GANDOLF_LIBRARY}")
# message("
# ")

##---------------------------------------------------------------------------##
## End of src/CMakeLists.txt
##---------------------------------------------------------------------------##
