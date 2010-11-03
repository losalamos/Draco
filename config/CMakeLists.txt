#-----------------------------*-cmake-*----------------------------------------#
# file   config/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 Nov 2
# brief  Install the config directory in the target directory 
# note   © Copyright 2010, Los Alamos National Security, All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

cmake_minimum_required( VERSION 2.8 FATAL_ERROR )
project( config )

file( GLOB CMake_src *.cmake )
file( GLOB CMake_in  *.in    )

install( 
   FILES  ${CMake_src}  ${CMake_in}
   DESTINATION config )

#include( uninstall )
#provide_uninstall_target()


##---------------------------------------------------------------------------##
## End of config/CMakeLists.txt
##---------------------------------------------------------------------------##
