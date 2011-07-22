project(draco_doc_releases NONE)
cmake_minimum_required(VERSION 2.8)

# the documents
file( GLOB latex_documents draco-*.tex )
list( REMOVE_ITEM latex_documents
   ${PROJECT_SOURCE_DIR}/draco-6_2_0.tex )

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
