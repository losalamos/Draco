#-----------------------------*-cmake-*----------------------------------------#
# file   doc/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2011 July 20
# brief  Instructions for building doc level Makefile.
# note   © Copyright 2011 Los Alamos National Security, All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

cmake_minimum_required( VERSION 2.8 FATAL_ERROR )
project( draco_doc_dir )

# Set default compile environment:
find_package(LATEX)

# Helper macro that only tries to configure a listed package if it
# exists in the source tree.  For example, when building milagro, the
# get_draco script only checks out the parts of Draco that are
# required by milagro.
macro( add_dir_if_exists package )
   message(STATUS "${package}")
   if( EXISTS ${PROJECT_SOURCE_DIR}/${package} )
      add_subdirectory( ${package} )
   endif()
endmacro()

# Add directories
add_subdirectory( releases )

##---------------------------------------------------------------------------##
## End of doc/CMakeLists.txt
##---------------------------------------------------------------------------##
