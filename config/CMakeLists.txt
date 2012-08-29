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

file( GLOB CMake_src *.cmake )
file( GLOB CMake_in  *.in    )

# include(debug_macros)
# echo_all_cmake_variable_values()

string( TOUPPER "${CMAKE_BUILD_TYPE}" upper_build_type )
set( Draco_CXX_COMPILER_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${upper_build_type}}" )
set( Draco_C_COMPILER_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${upper_build_type}}" )
set( Draco_Fortran_COMPILER_FLAGS "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${upper_build_type}}" )

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
