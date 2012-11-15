#-----------------------------*-cmake-*----------------------------------------#
# file   config/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 Nov 2
# brief  Install the config directory in the target directory 
# note   Copyright (C) 2010-2012, Los Alamos National Security
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

cmake_minimum_required( VERSION 2.8 FATAL_ERROR )
project( config )

set( CMake_src
#   FindGSL.cmake
#   FindGrace.cmake
#   FindLIBSCI.cmake
#   FindRandom123.cmake
#   Toolchain-catamount.cmake
   Toolchain-roadrunner-ppe.cmake
   Toolchain-roadrunner-spu.cmake
   apple-clang.cmake
   buildEnv.cmake
   compilerEnv.cmake
   component_macros.cmake
   configureFileOnMake.cmake
   debug_macros.cmake
   doc_macros.cmake
   dracoTesting.cmake
   dracoVersion.cmake
   parse_arguments.cmake
   platform_checks.cmake
   platform_customization_roadrunner_ppe.cmake
   unix-g++.cmake
   unix-gfortran.cmake
   unix-ifort.cmake
   unix-intel.cmake
   unix-pgf90.cmake
   unix-pgi.cmake
   unix-ppu.cmake
   unix-spu.cmake
   unix-xl.cmake
   unix-xlf.cmake
#   vendor_libraries.cmake
   windows-cl.cmake
   windows-ifort.cmake
)
set( CMake_in
   cmake_uninstall.cmake.in
   doxygen_config.in
   )

# include(debug_macros)
# echo_all_cmake_variable_values()

string( TOUPPER "${CMAKE_BUILD_TYPE}" upper_build_type )
set( Draco_CXX_COMPILER_FLAGS 
   "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${upper_build_type}}" )
set( Draco_C_COMPILER_FLAGS 
   "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${upper_build_type}}" )
set( Draco_Fortran_COMPILER_FLAGS 
   "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${upper_build_type}}" )

# sanity checks
# 'src' directory must be processed before the 'config' directory.
if( "${Draco_PACKAGE_LIST}" STREQUAL "" )
   message( FATAL_ERROR 
"Something is wrong with the build system.  Draco_PACKAGE_LIST is empty." )
endif()
if( "${Draco_TPL_LIST}" STREQUAL "" )
   message( FATAL_ERROR 
"Something is wrong with the build system.  Draco_TPL_LIST is empty." )
endif()
if( "${Draco_LIBRARIES}" STREQUAL "" )
   message( FATAL_ERROR 
"Something is wrong with the build system.  Draco_LIBRARIES is empty." )
endif()
if( "${Draco_TPL_LIBRARIES}" STREQUAL "" )
   message( FATAL_ERROR 
"Something is wrong with the build system.  Draco_TPL_LIBRARIES is empty." )
endif()

# Process the config file
configure_file( draco-config-install.cmake.in
   ${Draco_BINARY_DIR}/CMakeFiles/draco-config.cmake @ONLY)

install( 
   FILES  
     ${CMake_src}  
     ${CMake_in} 
     ${Draco_BINARY_DIR}/CMakeFiles/draco-config.cmake
   DESTINATION
     lib/cmake/draco-${DRACO_VERSION_MAJOR}.${DRACO_VERSION_MINOR} )

##---------------------------------------------------------------------------##
## End of config/CMakeLists.txt
##---------------------------------------------------------------------------##
