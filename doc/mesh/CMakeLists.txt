#-----------------------------*-cmake-*----------------------------------------#
# file   draco/doc/matprops/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 October 29
# brief  Instructions for building LaTeX documentation.
# note   Copyright (C) 2012 Los Alamos National Security, 
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id: CMakeLists.txt 6817 2012-10-11 18:08:03Z kellyt $
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 2.8)
project(dracodoc_mesh NONE)

if( "${PROJECT_SOURCE_DIR}" MATCHES "${PROJECT_BINARY_DIR}" )
   message( FATAL_ERROR "Must use out-of-source build tree!" )
endif()

# Build system configuration files are located here.
if( NOT "${CMAKE_MODULE_PATH}" MATCHES "${PROJECT_SOURCE_DIR}/../../config" )
   set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../../config )
endif()
include( UseLATEX )

##---------------------------------------------------------------------------##

add_latex_document( Expression.tex
   DEFAULT_PDF
   BIBFILES ../../environment/bibfiles/draco.bib
   MANGLE_TARGET_NAMES
)

add_latex_document( MTConcept.tex
   DEFAULT_PDF
   BIBFILES ../../environment/bibfiles/draco.bib
   MANGLE_TARGET_NAMES
)