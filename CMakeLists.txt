#-----------------------------*-cmake-*----------------------------------------#
# file   draco/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 April 28
# brief  Instructions for building root level Makefile.
# note   Copyright (C) 2010-2013 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# How to build Draco using CMake:
#------------------------------------------------------------------------------#
# Background:
# 
# There are 3 directories: source; build; target.
#   source - this is the location of svn co files
#   build  - this is the location for makefiles, object files, etc.
#   target - this is the install location for final headers, libraries
#            and binaries. 
#
# Environment:
#
# - CMake must be in your path.
# - The compiler is chosen by setting $CXX, $CC, $FC, etc.  On the
#   ICN this is done for you automatically.  If not set, CMake will
#   use whatever it can find in your path (e.g.: /usr/bin/g++)
# - Vendor locations can be specified by setting $<VENDOR>_LIB_DIR and
#   $<VENDOR>_INC_DIR. For example, $MPI_INC_DIR, $GSL_LIB_DIR,
#   etc. When using the Draco modules feature, this is done
#   automatically. 
#
# Requirements:
#
# 1. By default MPI is selected by finding 
#    mpirun/mpiexec/aprun/srun/etc in your PATH. 
# 2. GSL must be found in a standard system location or the
#    environment variables GSL_INC_DIR and GSL_LIB_DIR must be set to
#    point to the local GSL installation.
#
# Building:
# 
# 1. Change to your build directory.
# 2. [optional] Copy the $SRC_DIR/CMakeCache.txt file to the build
#    directory and edit the file as needed.
# 3. Run 'cmake $SRC_DIR'.  This generates Makefiles in the build
#    directory. 
# 4. [optional] Use command line arguments to cmake to configure your
#    build (instead of using CMakeCache.txt).  For example:
#    cmake -DCMAKE_INSTALL_PREFIX=/var/tmp/draco-install-dir \
#          -DDRACO_DBC_LEVEL=0 \
#          -DCMAKE_BUILD_TYPE=Release \
#          /path/to/draco/source
# 5. [optional] Run 'make edit_cache' to review your build setup.
#    Edit as needed. 
# 6. Run 'make -j' to compile the code. (also try 'make VERBOSE=1')
# 7. Run 'ctest -j16' to run the tests.
# 8. Run 'make install' to install libraries in the $CMAKE_PREFIX_DIR.
# 9. [optional] Run 'make help' for a list of targets.
#
# GUI
#
# If you prefer working in a GUI instead of using the command line,
# run 'cmake-gui' and provide values for 'where is the source code'
# and 'where to build the binaries'.  Press the 'configure' and
# 'generate' buttons to generate Makefiles or project files. 
#
# Development Environments
#
# If you prefer working in XCode or Eclipse CDT, you should tell CMake
# to generate project files instead of makefiles:
#    cmake -G"Eclipse CDT4 - Unix Makefiles" ...
#    cmake -G"XCode"
#
# Reference:
#
# Draco build system: /ccs/codes/radtran/svn/draco/dracodoc/build or
#    https://rtt.lanl.gov/redmine/projects/draco/repository/dracodoc/changes/build/dbs.pdf
# CMake Tutorial: http://www.cmake.org/cmake/help/cmake_tutorial.html
# CMake Wiki:     http://www.cmake.org/Wiki/CMake
# CMake FAQ:      http://www.cmake.org/Wiki/CMake_FAQ
# CMake Documentation: http://www.cmake.org/cmake/help/documentation.html
# CMake version compatibility matrix: 
#       http://www.cmake.org/Wiki/CMake_Version_Compatibility_Matrix
#------------------------------------------------------------------------------#

cmake_minimum_required( VERSION 2.8.8 FATAL_ERROR ) 
# 2.8.7 Improved support for Eclipse CDT.  Add support for CUDA_PATH.
# 2.8.8 Support for CMAKE_CXX_COMPILER_VERSION.
#       Addition of CMakePackageConfigHelpers
#       Support MSVC builds mixed with MingGW gfortran
#       Support for Xcode 4.3
project( Draco CXX C )
 
 # Do not look for Fortran for 
 # 1. XCode based Generators, or
 # 2. Visual Studio IDE or NMake Generators (MSYS or CYGWIN environments will
 #    look for Fortran. 
 # 3. Ninja, Codeblocks or Eclipse CDT4 generators.
if( ${CMAKE_GENERATOR} MATCHES "Unix Makefiles" )
   enable_language( Fortran OPTIONAL )
endif()


# Build system configuration files are located here.
set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config )

#
# Debug cmake scripts:
#
# 1. Review macros defined at config/debug_macros.cmake.
# 2. Uncomment this feature to tell CMake to print include paths
#    during target registration (requires cmake-2.8.11+).
# set( CMAKE_DEBUG_TARGET_PROPERTIES INCLUDE_DIRECTORIES )

#
# The Draco version number.
#
set(DRACO_VERSION_MAJOR 6)
set(DRACO_VERSION_MINOR 10)
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

# Set compiler options
include( compilerEnv )
dbsSetupCompilers()
dbsSetupCxx()

# Is Fortran enabled (it is considered 'optional' for draco)?
get_property(_LANGUAGES_ GLOBAL PROPERTY ENABLED_LANGUAGES)
if( _LANGUAGES_ MATCHES Fortran )
   dbsSetupFortran()
endif()

include( FeatureSummary )

# Find any globally required libraries
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
   FILES 
      ChangeLog
      Copyright
      README.draco
      ToDo
   DESTINATION ${CMAKE_INSTALL_PREFIX}
)

# Add support for the 'uninstall' target
configure_file(
   "${CMAKE_CURRENT_SOURCE_DIR}/config/cmake_uninstall.cmake.in"
   "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
   @ONLY )
add_custom_target( uninstall
   "${CMAKE_COMMAND}" -P
   "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake" )
# Sometimes CMake does not install the correct RPATH.  This target
# usually fixes the issue
add_custom_target( update_rpaths
  COMMAND "${CMAKE_MAKE_PROGRAM}" rebuild_cache
  COMMAND "${CMAKE_MAKE_PROGRAM}" ${MAKE_ARGS} install
  COMMENT "Reparsing to try to update RPATH values..."
  )

# Export targets
install( EXPORT draco-targets DESTINATION lib/cmake/draco )

##---------------------------------------------------------------------------##
## Notes:

# 1. Use 'cmake -LAH' to list all user set variables.
# 2. Use 'make rebuild_cache' to regenerate config.h files, directories,
#   etc.
# 3. Use 'make edit_cache' to change user set variables.
# 4. Use 'make -jN' to build in parallel with N processes.
# 5. Use 'make test' or 'ctest -jN' to run the unit tests.
# 6. Use 'make install' to strip libraries and copy libraries and
#    headers to the INSTALL_PREFIX directory.

##---------------------------------------------------------------------------##
## End of CMakeLists.txt
##---------------------------------------------------------------------------##
