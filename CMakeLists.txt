#-----------------------------*-cmake-*----------------------------------------#
# file   draco/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 April 28
# brief  Instructions for building root level Makefile.
# note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

# Build notes: http://rtt/redmine/projects/draco/wiki/Configuration_and_Building

cmake_minimum_required(VERSION 3.5.0 FATAL_ERROR)
project( Draco CXX C )

# Do not look for Fortran for
# 1. XCode based Generators, or
# 2. Visual Studio IDE or NMake Generators (MSYS or CYGWIN environments will
#    look for Fortran).
# 3. Ninja, Codeblocks or Eclipse CDT4 generators.
# 4. Unix Makefile generators when $ENV{FC} is not set (clang on Linux)
if( ${CMAKE_GENERATOR} MATCHES "Unix Makefiles" AND
    ( NOT "$ENV{FC}x" STREQUAL "x" OR CMAKE_Fortran_COMPILER ))
  enable_language( Fortran OPTIONAL )
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
set(DRACO_VERSION_MAJOR 6)
set(DRACO_VERSION_MINOR 22)
# set(DRACO_VERSION_PATCH 0) # only set this for a release
include(dracoVersion)
set_ccs2_software_version( DRACO )

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
message("
Platform Checks...
")
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
dbsSetupProfilerTools()

# Find any globally required libraries
include( FeatureSummary )
include( vendor_libraries )
setupVendorLibraries()

#
# Build Draco components:
#
# add_subdirectory( doc ) Now at svnroot/draco/dracodoc.
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
  FILES ChangeLog Copyright README.draco
  DESTINATION ${CMAKE_INSTALL_PREFIX} )

# Add support for the 'uninstall' target
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/config/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  @ONLY )
add_custom_target( uninstall
  COMMAND "${CMAKE_COMMAND}" -P
          "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake" )
# Sometimes CMake does not install the correct RPATH.  This target
# usually fixes the issue
add_custom_target( update_rpaths
  COMMAND "${CMAKE_MAKE_PROGRAM}" rebuild_cache
  COMMAND "${CMAKE_MAKE_PROGRAM}" ${MAKE_ARGS} install
  COMMENT "Reparsing to try to update RPATH values..."
  )

# Export targets
install(
  EXPORT draco-targets
  DESTINATION cmake/draco
  EXPORT_LINK_INTERFACE_LIBRARIES )

##---------------------------------------------------------------------------##
## End of CMakeLists.txt
##---------------------------------------------------------------------------##
