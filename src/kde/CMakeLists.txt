#--------------------------------------------*-cmake-*---------------------------------------------#
#
# file   kde/CMakeLists.txt
# author Mathew Cleveland
# date   Nov. 10th 2020
# brief  CMake build environment for KDE
# note   Copyright (C) 2018-2020 Triad National Security, LLC.
#        All rights reserved.
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.9.0)
project( kde CXX )

#--------------------------------------------------------------------------------------------------#
# Source files
#--------------------------------------------------------------------------------------------------#

file( GLOB template_implementations *.t.hh *.i.hh )
file( GLOB sources *.cc )
file( GLOB explicit_instantiations *_pt.cc )
file( GLOB headers *.hh )
list( REMOVE_ITEM headers ${template_implementations} )

# Make the header files available in the IDE.
if( MSVC_IDE OR ${CMAKE_GENERATOR} MATCHES Xcode )
   list( APPEND sources ${headers} ${template_implementations} )
endif()

#--------------------------------------------------------------------------------------------------#
# Directories to search for include directives
#--------------------------------------------------------------------------------------------------#

include_directories(
   ${PROJECT_SOURCE_DIR}   # sources
   ${PROJECT_BINARY_DIR} ) # config.h

#--------------------------------------------------------------------------------------------------#
# Build package library
#--------------------------------------------------------------------------------------------------#

add_component_library(
   TARGET       Lib_kde
   TARGET_DEPS  Lib_dsxx Lib_c4
   LIBRARY_NAME ${PROJECT_NAME}
   HEADERS      "${headers}"
   SOURCES      "${sources}"
   )

#--------------------------------------------------------------------------------------------------#
# Installation instructions
#--------------------------------------------------------------------------------------------------#

install( FILES ${headers} ${template_implementations}
  DESTINATION ${DBSCFGDIR}include/kde )

#--------------------------------------------------------------------------------------------------#
# Unit tests
#--------------------------------------------------------------------------------------------------#

if( BUILD_TESTING )
   add_subdirectory( test )
endif()

#--------------------------------------------------------------------------------------------------#
# Autodoc
#--------------------------------------------------------------------------------------------------#

process_autodoc_pages()
#--------------------------------------------------------------------------------------------------#
# End kde/CMakeLists.txt
#--------------------------------------------------------------------------------------------------#

