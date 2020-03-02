#-----------------------------*-cmake-*----------------------------------------#
# file   config/doc_macros.cmake
# author Kelly Thompson
# date   2011 July 22
# brief  Provide extra macros to simplify CMakeLists.txt for latex doc
#        directories.
# note   Copyright (C) 2016-2020 Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

# setup latex environment
find_package(LATEX)

##---------------------------------------------------------------------------##
## Build pdf from latex
#
# Given a list of latex source master files create pdf documents, one
# for each listed latex source file.  Create a makefile target that
# will be ${docname}.pdf
#
# Usage:
#
# add_latex_doc(
#    SOURCES   "${list_of_latex_files}"
#    TEXINPUTS "${PROJECT_SOURCE_DIR}:${PROJECT_SOURCE_DIR}/../../environment/latex"
#    BIBINPUTS "${PROJECT_SOURCE_DIR}/../../environment/bibfiles"
# )
#
##---------------------------------------------------------------------------##
macro( add_latex_doc )

   # These become variables of the form ${addscalartests_SOURCES}, etc.
   cmake_parse_arguments(
      # prefix
      ald
      # list names
      "SOURCES;BIBINPUTS;TEXINPUTS"
      # option names
      "NONE"
      ${ARGV}
      )

   set( ald_BIBINPUTS "${ald_BIBINPUTS}:$ENV{BSTINPUTS}" )
   set( ald_TEXINPUTS "${ald_TEXINPUTS}:$ENV{TEXINPUTS}" )

   # compile commands
   if(LATEX_COMPILER)

      # if LATEX_COMPILER assume that DVIPS_CONVERTER and
      # PS2PDF_CONVERTER are also available.

      foreach( fpdoc ${ald_SOURCES} )
         get_filename_component( doc ${fpdoc} NAME_WE )
         add_custom_command(
            OUTPUT    ${PROJECT_BINARY_DIR}/${doc}.dvi
            COMMAND   TEXINPUTS=${ald_TEXINPUTS}
                      ${LATEX_COMPILER}  ${PROJECT_SOURCE_DIR}/${doc}.tex
            COMMAND   BIBINPUTS=${ald_BIBINPUTS}
                      ${BIBTEX_COMPILER} ${PROJECT_BINARY_DIR}/${doc}.aux
            COMMAND   TEXINPUTS=${ald_TEXINPUTS}
                      ${LATEX_COMPILER} ${PROJECT_SOURCE_DIR}/${doc}.tex
            COMMAND   TEXINPUTS=${ald_TEXINPUTS}
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
         add_custom_target(${doc}.pdf # ALL
            DEPENDS   ${PROJECT_BINARY_DIR}/${doc}.pdf
            )

         # target to build all pdf files
         list( APPEND dirpdfs ${PROJECT_BINARY_DIR}/${doc}.pdf )
         list( APPEND dirpdfs_tgt ${doc}.pdf )

         # generate a list of extra files to remove during 'make clean'
         list( APPEND extra_files
            ${PROJECT_BINARY_DIR}/${doc}.log
            ${PROJECT_BINARY_DIR}/${doc}.bbl
            ${PROJECT_BINARY_DIR}/${doc}.blg
            ${PROJECT_BINARY_DIR}/${doc}.aux )

      endforeach()

      # Provide extra instructions for 'make clean'
      set_directory_properties(
         PROPERTIES
         ADDITIONAL_MAKE_CLEAN_FILES "${extra_files}" )

      # Create a target to generate all pdfs for this directory.
      add_custom_target(${PROJECT_NAME}_pdf # ALL
         DEPENDS ${dirpdfs} )

      get_directory_property( pdir PARENT_DIRECTORY )
      if( EXISTS "${pdir}" )
         set( allpdfs ${allpdfs};${dirpdfs_tgt} PARENT_SCOPE )
      endif()

      # Install instructions
      install( FILES ${dirpdfs} DESTINATION doc OPTIONAL)

   endif(LATEX_COMPILER)

endmacro()


##---------------------------------------------------------------------------##
## Build pdf from latex
#
# Given a list of latex source master files create pdf documents, one
# for each listed latex source file.  Create a makefile target that
# will be ${docname}.pdf
#
# Usage:
#
# add_latex_doc(
#    SOURCES   "${list_of_latex_files}"
#    TEXINPUTS "${PROJECT_SOURCE_DIR}:${PROJECT_SOURCE_DIR}/../../environment/latex"
#    BIBINPUTS "${PROJECT_SOURCE_DIR}/../../environment/bibfiles"
# )
#
##---------------------------------------------------------------------------##
macro( add_pdflatex_doc )

   # These become variables of the form ${addscalartests_SOURCES}, etc.
   cmake_parse_arguments(
      # prefix
      ald
      # list names
      "SOURCES;BIBINPUTS;TEXINPUTS"
      # option names
      "NONE"
      ${ARGV}
      )

   set( ald_BIBINPUTS "${ald_BIBINPUTS}:$ENV{BSTINPUTS}" )
   set( ald_TEXINPUTS "${ald_TEXINPUTS}:$ENV{TEXINPUTS}" )

   # compile commands
   if(PDFLATEX_COMPILER)

      foreach( fpdoc ${ald_SOURCES} )
         get_filename_component( doc ${fpdoc} NAME_WE )
         add_custom_command(
            OUTPUT    ${PROJECT_BINARY_DIR}/${doc}.pdf
            COMMAND   TEXINPUTS=${ald_TEXINPUTS}
                      ${PDFLATEX_COMPILER} ${PROJECT_SOURCE_DIR}/${doc}.tex
            COMMAND   BIBINPUTS=${ald_BIBINPUTS}
                      ${BIBTEX_COMPILER} ${PROJECT_BINARY_DIR}/${doc}.aux
            COMMAND   TEXINPUTS=${ald_TEXINPUTS}
                      ${PDFLATEX_COMPILER} ${PROJECT_SOURCE_DIR}/${doc}.tex
            COMMAND   TEXINPUTS=${ald_TEXINPUTS}
                      ${PDFLATEX_COMPILER} ${PROJECT_SOURCE_DIR}/${doc}.tex
            DEPENDS   ${PROJECT_SOURCE_DIR}/${doc}.tex
            COMMENT   "pdflatex (${doc}.tex --> ${doc}.pdf)"
            )

         # register target for each individual pdf file (eg: make draco-3_0_0.pdf)
         add_custom_target(${doc}.pdf # ALL
            DEPENDS   ${PROJECT_BINARY_DIR}/${doc}.pdf
            )

         # target to build all pdf files
         list( APPEND dirpdfs ${PROJECT_BINARY_DIR}/${doc}.pdf )
         list( APPEND dirpdfs_tgt ${doc}.pdf )

         # generate a list of extra files to remove during 'make clean'
         list( APPEND extra_files
            ${PROJECT_BINARY_DIR}/${doc}.log
            ${PROJECT_BINARY_DIR}/${doc}.bbl
            ${PROJECT_BINARY_DIR}/${doc}.blg
            ${PROJECT_BINARY_DIR}/${doc}.aux )

      endforeach()

      # Provide extra instructions for 'make clean'
      set_directory_properties(
         PROPERTIES
         ADDITIONAL_MAKE_CLEAN_FILES "${extra_files}" )

      # Create a target to generate all pdfs for this directory.
      add_custom_target(${PROJECT_NAME}_pdflatex # ALL
         DEPENDS ${dirpdfs} )

      get_directory_property( pdir PARENT_DIRECTORY )
      if( EXISTS "${pdir}" )
         set( allpdfs ${allpdfs};${dirpdfs_tgt} PARENT_SCOPE )
      endif()

      # Install instructions
      install( FILES ${dirpdfs} DESTINATION doc OPTIONAL )

   endif(PDFLATEX_COMPILER)

endmacro()
