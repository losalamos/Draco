#-----------------------------*-cmake-*----------------------------------------#
# file   config/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2013 June 24
# brief  Generate build project files for fortran.
# note   Copyright (C) 2013-2014, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id: CMakeLists.txt 6705 2012-08-29 16:21:06Z kellyt $
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 2.8)
project( FortranChecks_f90sub Fortran )

# ---------------------------------------------------------------------------- #
# Setup - This folder is parsed independent of parent folders, so we need 
# some basic setup.
# ---------------------------------------------------------------------------- #
# Build system configuration files are located here.
if( NOT EXISTS ${PROJECT_SOURCE_DIR}/../../../config )
   message( FATAL_ERROR "can't find ${PROJECT_SOURCE_DIR}/../../../config" )
else()
   set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../../../config )
endif()

# Provide helper functions used by component CMakeLists.txt files.
include( component_macros )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #
file( GLOB f90sources *.f90 )

# ---------------------------------------------------------------------------- #
# Build package libraries
# ---------------------------------------------------------------------------- #
add_component_library(
   TARGET       Lib_FC_f90sub
   LIBRARY_NAME FC_f90sub
   SOURCES      "${f90sources}" 
   LINK_LANGUAGE "Fortran" )
