project(draco_doc_releases NONE)
cmake_minimum_required(VERSION 2.8)

# Configuration
if( NOT "${CMAKE_MODULE_PATH}" MATCHES "config" )
   set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../../config )
   include( doc_macros )
   set( allpdfs "" CACHE INTERNAL "list of all pdfs to be generated" )
endif()

# the documents
file( GLOB latex_documents draco-*.tex )
set( pdflatex_documents 
  ${PROJECT_SOURCE_DIR}/draco-6_0_0.tex
  ${PROJECT_SOURCE_DIR}/draco-6_2_0.tex )
list( REMOVE_ITEM latex_documents ${pdflatex_documents} )

# LaTeX supporting files
# draco.bib is at draco/environment/bibfiles
# extra inputs (eps, etc) are at the source directory.
# cls files are at draco/environment/latex
set( BIBINPUTS "${PROJECT_SOURCE_DIR}/../../environment/bibfiles" )
set( TEXINPUTS
   "${PROJECT_SOURCE_DIR}:${PROJECT_SOURCE_DIR}/../../environment/latex"
   )

# Register LaTeX documents

add_latex_doc(
   SOURCES   "${latex_documents}" 
   BIBINPUTS "${BIBINPUTS}"
   TEXINPUTS "${TEXINPUTS}"
)

add_pdflatex_doc(
   SOURCES   "${pdflatex_documents}"
   BIBINPUTS "${BIBINPUTS}"
   TEXINPUTS "${TEXINPUTS}"
)

# Create a default target for this directory
if( NOT EXISTS ${draco_doc_dir_SOURCE_DIR} )
   add_custom_target( pdf ALL 
      DEPENDS ${PROJECT_BINARY_DIR}/draco-6_2_0.pdf )
endif()