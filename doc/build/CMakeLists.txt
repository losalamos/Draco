#-----------------------------*-cmake-*----------------------------------------#
# file   draco/doc/build/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 April 28
# brief  Instructions for building LaTeX documentation.
# note   Copyright (C) 2012 Los Alamos National Security, 
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id: CMakeLists.txt 6817 2012-10-11 18:08:03Z kellyt $
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 2.8)
project(dracodoc_buildsystem NONE)

if( "${PROJECT_SOURCE_DIR}" MATCHES "${PROJECT_BINARY_DIR}" )
   message( FATAL_ERROR "Must use out-of-source build tree!" )
endif()

# Build system configuration files are located here.
if( NOT "${CMAKE_MODULE_PATH}" MATCHES "${PROJECT_SOURCE_DIR}/../../config" )
   set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../../config )
endif()
include( UseLATEX )

##---------------------------------------------------------------------------##
set( main_latex_file dbs.tex )
set( image_files
   fig/build_tree.pdf
   fig/build_tree_post_install.pdf
   fig/cmakegui_defaults.pdf
   fig/level.pdf
   fig/package.pdf
   fig/src_config.pdf
   fig/src_doc.pdf
   fig/src_draco.pdf
   fig/src_src.pdf
   fig/src_tools.pdf
   fig/src_tree.pdf
)
set( extra_input_files
   adding.tex
   bashprofile.tex
   bashrc.tex
   compile.tex
   dbs.tex
   extend.tex
   extern.tex
   intro.tex
   model.tex
   quickstart.tex
   vendors.tex
)

add_latex_document( ${main_latex_file}
   INPUTS ${extra_input_files}
   IMAGES ${image_files}
   DEFAULT_PDF
   BIBFILES ../../environment/bibfiles/draco.bib
   USE_INDEX
)