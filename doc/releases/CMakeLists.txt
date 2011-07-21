project(draco_doc_releases NONE)
cmake_minimum_required(VERSION 2.8)

# the documents
file( GLOB latex_documents draco-*.tex )

# LaTeX supporting files
# draco.bib is at draco/environment/bibfiles
# extra inputs (eps, etc) are at the source directory.
# cls files are at draco/environment/latex
set( BSTINPUTS "${PROJECT_SOURCE_DIR}/../../environment/bibfiles:$ENV{BSTINPUTS}" )
set( TEXINPUTS "${PROJECT_SOURCE_DIR}:${PROJECT_SOURCE_DIR}/../../environment/latex:$ENV{TEXINPUTS}" )

# setup latex environment
find_package(LATEX)

# compile commands
if(LATEX_COMPILER)

   # if LATEX_COMPILER assume that DVIPS_CONVERTER and
   # PS2PDF_CONVERTER are also available.

   foreach( fpdoc ${latex_documents} )
      get_filename_component( doc ${fpdoc} NAME_WE )
      add_custom_command( 
         OUTPUT    ${PROJECT_BINARY_DIR}/${doc}.dvi
         COMMAND   TEXINPUTS=${TEXINPUTS}
                   ${LATEX_COMPILER}  ${PROJECT_SOURCE_DIR}/${doc}.tex
         COMMAND   BIBINPUTS=${BSTINPUTS}
                   ${BIBTEX_COMPILER} ${PROJECT_BINARY_DIR}/${doc}.aux
         COMMAND   TEXINPUTS=${TEXINPUTS}
                   ${LATEX_COMPILER} ${PROJECT_SOURCE_DIR}/${doc}.tex
         COMMAND   TEXINPUTS=${TEXINPUTS}
                   ${LATEX_COMPILER} ${PROJECT_SOURCE_DIR}/${doc}.tex
         DEPENDS   ${PROJECT_SOURCE_DIR}/${doc}.tex
         COMMENT   "Tex2dvi (${doc}.tex --> ${doc}.dvi)"
         )
      
      add_custom_command( 
         OUTPUT    ${PROJECT_BINARY_DIR}/${doc}.ps
         COMMAND   ${DVIPS_CONVERTER}
                   ${PROJECT_BINARY_DIR}/${doc}.dvi
                   -o ${PROJECT_BINARY_DIR}/${doc}.ps
         DEPENDS   ${PROJECT_BINARY_DIR}/${doc}.dvi
         COMMENT   "dvi2ps (${doc}.dvi --> ${doc}.ps)"
         )
      
      
      add_custom_command( 
         OUTPUT    ${PROJECT_BINARY_DIR}/${doc}.pdf
         COMMAND   ${PS2PDF_CONVERTER}
                   ${PROJECT_BINARY_DIR}/${doc}.ps
         DEPENDS   ${PROJECT_BINARY_DIR}/${doc}.ps
         COMMENT   "ps2pdf (${doc}.ps --> ${doc}.pdf)"
         )
      
      # register target for each individual pdf file (eg: make draco-3_0_0.pdf)
      add_custom_target(${doc}.pdf ALL 
         DEPENDS   ${PROJECT_BINARY_DIR}/${doc}.pdf
         )
      
      # target to build all pdf files
      list( APPEND allpdf
         ${PROJECT_BINARY_DIR}/${doc}.pdf )
      
      # generate a list of extra files to remove during 'make clean'
      list( APPEND extra_files 
         ${PROJECT_BINARY_DIR}/${doc}.log
         ${PROJECT_BINARY_DIR}/${doc}.bbl
         ${PROJECT_BINARY_DIR}/${doc}.blg
         ${PROJECT_BINARY_DIR}/${doc}.aux )
      
   endforeach()

   # Create a target to generate all pdfs for this directory.
   add_custom_target(${PROJECT_NAME}_pdf ALL DEPENDS ${allpdf} )
   
endif(LATEX_COMPILER)

# Provide extra instructions for 'make clean'
set_directory_properties(
   PROPERTIES
   ADDITIONAL_MAKE_CLEAN_FILES "${extra_files}" )

# Install instructions
install( FILES ${allpdf} DESTINATION doc )
