#--------------------------------------------*-cmake-*---------------------------------------------#
# file   draco/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 April 28
# brief  Instructions for building root level Makefile.
# note   Copyright (C) 2019-2020 Triad National Security, LLC.
#        All rights reserved.
#--------------------------------------------------------------------------------------------------#

# Build notes:
#   https://rtt.lanl.gov/redmine/projects/draco/wiki/
#   https://rtt.lanl.gov/redmine/projects/draco/wiki/Common_Configure_Options

cmake_minimum_required(VERSION 3.17.0 FATAL_ERROR)
set(ddesc "An object-oriented component library supporting radiation")
string(APPEND ddesc " transport applications.")
project( Draco DESCRIPTION ${ddesc} VERSION 7.8 LANGUAGES CXX C)

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
  if( DEFINED ENV{CUDADIR} OR DEFINED ENV{CUDACXX} OR DEFINED ENV{CUDA_HOME} )
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
dbsSetupCuda()
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
if( TARGET Exe_draco_info )
  add_subdirectory( autodoc ) # This must be processed after 'src'
endif()
add_subdirectory( config )  # This must be processed after 'src'

# install top level documents
install(
  FILES ChangeLog LICENSE.md README.md
  DESTINATION ${CMAKE_INSTALL_PREFIX} )

# Export targets
install(
  EXPORT draco-targets
  DESTINATION cmake
  EXPORT_LINK_INTERFACE_LIBRARIES )

##---------------------------------------------------------------------------##
## End of CMakeLists.txt
##---------------------------------------------------------------------------##
