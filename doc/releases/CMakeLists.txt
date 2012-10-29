#-----------------------------*-cmake-*----------------------------------------#
# file   draco/doc/releases/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 October 29
# brief  Instructions for building LaTeX documentation.
# note   Copyright (C) 2012 Los Alamos National Security, 
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 2.8)
project(dracodoc_releases NONE)

if( "${PROJECT_SOURCE_DIR}" MATCHES "${PROJECT_BINARY_DIR}" )
   message( FATAL_ERROR "Must use out-of-source build tree!" )
endif()

# Build system configuration files are located here.
if( NOT "${CMAKE_MODULE_PATH}" MATCHES "${PROJECT_SOURCE_DIR}/../../config" )
   set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../../config )
endif()
include( UseLATEX )

##---------------------------------------------------------------------------##

set( image_files
   dashboard-620.jpg
   dashboard-cc-620.jpg
   dashboard-valgrind-620.jpg
   loc-6_2_0.jpg
   loc-3_0_0.eps
   loc-3_2_0.eps
   loc-4_0_0.eps
   loc-5_0_0.eps
   draco600-loc-history.pdf
   level-3_0_0.pdf
   level-5_0_0.pdf
   level-6_0_0.pdf
   level-6_2_0.pdf
   loc-3_0_0.pdf
   loc-3_2_0.pdf
   loc-4_0_0.pdf
   loc-5_0_0.pdf
   )

set( extra_inputs
   cov-5_0_0.tex
   cov-6_0_0.tex
   cov-6_2_0.tex
)

add_latex_document( draco-6_2_0.tex
   INPUTS ${extra_inputs}
   DEFAULT_PDF
   BIBFILES ../../environment/bibfiles/draco.bib
   IMAGES ${image_files}
   MANGLE_TARGET_NAMES
)

add_latex_document( draco-6_0_0.tex
   INPUTS ${extra_inputs}
   DEFAULT_PDF
   BIBFILES ../../environment/bibfiles/draco.bib
   IMAGES ${image_files}
   MANGLE_TARGET_NAMES
)

add_latex_document( draco-5_0_0.tex
   INPUTS ${extra_inputs}
   DEFAULT_PDF
   BIBFILES ../../environment/bibfiles/draco.bib
   IMAGES ${image_files}
   MANGLE_TARGET_NAMES
)

add_latex_document( draco-4_0_0.tex
   INPUTS ${extra_inputs}
   DEFAULT_PDF
   BIBFILES ../../environment/bibfiles/draco.bib
   IMAGES ${image_files}
   MANGLE_TARGET_NAMES
)

add_latex_document( draco-3_2_0.tex
   INPUTS ${extra_inputs}
   DEFAULT_PDF
   BIBFILES ../../environment/bibfiles/draco.bib
   IMAGES ${image_files}
   MANGLE_TARGET_NAMES
)

add_latex_document( draco-3_0_0.tex
   INPUTS ${extra_inputs}
   DEFAULT_PDF
   BIBFILES ../../environment/bibfiles/draco.bib
   IMAGES ${image_files}
   MANGLE_TARGET_NAMES
)