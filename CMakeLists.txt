#-----------------------------*-cmake-*----------------------------------------#
# file   draco/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 April 28
# brief  Instructions for building root level Makefile.
# note   Copyright (C) 2016-2019 Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

# Build notes:
#   https://rtt.lanl.gov/redmine/projects/draco/wiki/
#   https://rtt.lanl.gov/redmine/projects/draco/wiki/Common_Configure_Options

cmake_minimum_required(VERSION 3.11.0 FATAL_ERROR)
project( Draco
  VERSION 7.1
  DESCRIPTION "An object-oriented component library supporting radiation transport applications."
  LANGUAGES CXX C )

# Do not look for Fortran/CUDA for
# 1. XCode based Generators, or
# 2. Visual Studio IDE or NMake Generators (MSYS or CYGWIN environments will
#    look for Fortran).
# 3. Ninja, Codeblocks or Eclipse CDT4 generators.
# 4. Unix Makefile generator when:
#    - Fortran only: $ENV{FC} is not set (e.g. clang on Linux)
#    - Cuda only   : $ENV{CUDADIR} is not set (e.g. cuda module not loaded)
if( ${CMAKE_GENERATOR} MATCHES "Unix Makefiles")
  if( DEFINED ENV{FC} OR DEFINED CMAKE_Fortran_COMPILER )
    enable_language( Fortran OPTIONAL )
  endif()
  if( DEFINED ENV{CUDADIR} )
    enable_language( CUDA OPTIONAL )
  endif()
endif()

# Build system configuration files are located here.
set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config )

#
# Debug cmake scripts:
#
# 1. Review macros defined at config/debug_macros.cmake and
#    config/print_target_properties.cmake.
# 2. Uncomment this feature to tell CMake to print include paths
#    during target registration (requires cmake-2.8.11+).
#    set( CMAKE_DEBUG_TARGET_PROPERTIES INCLUDE_DIRECTORIES )
# 3. Optionally use --trace

#
# The Draco version number.
#
include(dracoVersion)
set_ccs2_software_version( Draco )

#
# Unit Test Setup
#
include(dracoTesting)

# Set default compile environment:

# Setup defaults, value checks, etc.
include(buildEnv)
# set defaults for BUILD_TYPE and INSTALL_PREFIX
dbsSetDefaults()
# Initialize fields that define the exported target files
# (draco-config.cmake)
dbsInitExportTargets( "Draco" )
# Save config info
dbsConfigInfo()

# Platform Checks: Is HOST_NAME_MAX defined?  Is WinSock2.h available?  Is
# gethostname() available?
include( platform_checks )
query_craype()
set_draco_uname()
query_have_gethostname()
query_have_maxpathlen()
query_have_sys_headers() # sets HAVE_UNISTD_H, etc.
query_have_restrict_keyword()
query_fma_on_hardware()

# Set compiler options
include( compilerEnv )
dbsSetupCompilers()
dbsSetupCxx()
dbsSetupFortran()
dbsSetupProfilerTools()
dbsSetupStaticAnalyzers()

# Find any globally required libraries
include( FeatureSummary )
include( vendor_libraries )
setupVendorLibraries()

#
# Build Draco components:
#
add_subdirectory( src )
add_subdirectory( autodoc ) # This must be processed after 'src'
add_subdirectory( config )  # This must be processed after 'src'

#
# Dashboard/Regression setup
#
# CTest/CDash custom configuration file
configure_file(
  ${PROJECT_SOURCE_DIR}/CTestCustom.cmake
  ${PROJECT_BINARY_DIR}/CTestCustom.cmake
  COPYONLY )

# install top level documents
install(
  FILES ChangeLog LICENSE.md README.md
  DESTINATION ${CMAKE_INSTALL_PREFIX} )

# Add support for the 'uninstall' target
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/config/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  @ONLY )
add_custom_target( uninstall
  COMMAND "${CMAKE_COMMAND}" -P
          "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake" )

# Export targets
install(
  EXPORT draco-targets
  DESTINATION cmake/draco
  EXPORT_LINK_INTERFACE_LIBRARIES )

##---------------------------------------------------------------------------##
## End of CMakeLists.txt
##---------------------------------------------------------------------------##
