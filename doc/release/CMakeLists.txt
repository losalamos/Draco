# 
# How to build documentation:
#
# From the source directory draco/doc/release:
# % cmake .
# % make [all|clean|draco-release.pdf]
#-----------------------------------------------------------------------------#
project(draco_doc_releases NONE)
cmake_minimum_required(VERSION 2.8)

# Configuration
if( NOT "${CMAKE_MODULE_PATH}" MATCHES "config" )
   find_file( dracodocmacros
      doc_macros.cmake
      PATHS 
        ${Draco_SOURCE_DIR}/config
        ${PROJECT_SOURCE_DIR}/../../config
   )
   get_filename_component( draco_config_dir ${dracodocmacros} PATH )
   set( CMAKE_MODULE_PATH ${draco_config_dir} )
   include( doc_macros )
   set( allpdfs "" CACHE INTERNAL "list of all pdfs to be generated" )
endif()

# the documents
#file( GLOB latex_documents draco-*.tex )
file( GLOB pdflatex_documents draco-*.tex )
#set( pdflatex_documents ${PROJECT_SOURCE_DIR}/draco-release.tex )
#list( REMOVE_ITEM latex_documents ${pdflatex_documents} )

# LaTeX supporting files
# draco.bib is at draco/environment/bibfiles
# IMC.bib is at imcdoc/bib.  
# extra inputs (eps, etc) are at the source directory.
# cls files are at draco/environment/latex

# Assume imcdoc is a directory parallel to clubimc.
find_file( dracobib
   NAMES draco.bib
   PATHS ${Draco_SOURCE_DIR}/environment/bibfiles
         ${PROJECT_SOURCE_DIR}/../../environment/bibfiles
         $ENV{BIBINPUTS} 
         $ENV{BSTINPUTS}
         $ENV{TEXINPUTS}
         )
get_filename_component( draco_bib_dir ${dracobib} PATH )
if( NOT EXISTS ${draco_bib_dir} )
   message( FATAL_ERROR 
      "Expected directory draco/environment/bibfiles to be available.  This is required so that the bib file draco.bib can be found.")
endif()

set( BIBINPUTS "${DRACO_SOURCE_DIR}/environment/bibfiles:${draco_bib_dir}" )
set( TEXINPUTS
   "${PROJECT_SOURCE_DIR}:${DRACO_SOURCE_DIR}/environment/latex"
   )

# Register LaTeX documents

# add_latex_doc(
#    SOURCES   "${latex_documents}" 
#    BIBINPUTS "${BIBINPUTS}"
#    TEXINPUTS "${TEXINPUTS}"
# )

add_pdflatex_doc(
   SOURCES   "${pdflatex_documents}"
   BIBINPUTS "${BIBINPUTS}"
   TEXINPUTS "${TEXINPUTS}"
)

# default make target
add_custom_target( ${PROJECT_NAME}-pdf ALL
   DEPENDS ${PROJECT_BINARY_DIR}/draco-release.pdf )
