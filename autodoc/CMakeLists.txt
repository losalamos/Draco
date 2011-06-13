 #-----------------------------*-cmake-*----------------------------------------#
# file   draco/autodoc/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 Oct 14
# brief  Instructions for building Doxygen documentation
# note   Copyright © 2010 Los Alamos National Security, All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

# Doxygen

find_package( Doxygen QUIET )

if( DOXYGEN_FOUND )

  set( DOXYGEN_OUTPUT_DIR
    ${PROJECT_BINARY_DIR}/autodoc )

  # Source directories that should be examined by doxygen Sets
  # ${DOXYGEN_INPUT}.  This value is embedded in the doxygen_config
  # file below.
  set( DOXYGEN_INPUT 
     "${PROJECT_SOURCE_DIR}/autodoc ${PROJECT_BINARY_DIR}/autodoc" )

  file( GLOB package_list ${PROJECT_SOURCE_DIR}/src/* )
  foreach( package ${package_list} )
     if( EXISTS ${package}/CMakeLists.txt )
        set( DOXYGEN_INPUT "${DOXYGEN_INPUT} ${package}" )
        message( "Doxygen will process ${package}" )
     endif()
     if( EXISTS ${package}/test/CMakeLists.txt )
        set( DOXYGEN_EXAMPLES "${DOXYGEN_EXAMPLES} ${package}/test" )
        # message( "Doxygen will process examples found at ${package}/test" )
     endif()
     if( EXISTS ${package}/autodoc )
        set( DOXYGEN_INPUT "${DOXYGEN_INPUT} ${package}/autodoc" )
     endif()
  endforeach()  

  set( DOXYGEN_HTML_OUTPUT html )
  set( DOXYGEN_DOTFILE_DIRS dot )

  set( htmlfooter "html/footer.html")
  set( htmlheader "html/header.html")

# doxygen_config.in needs the following symbols:

  add_custom_command(
    OUTPUT  "${PROJECT_BINARY_DIR}/autodoc/doxygen_config"
    COMMAND "${CMAKE_COMMAND}" 
            -DINFILE="${PROJECT_SOURCE_DIR}/config/doxygen_config.in"
            -DOUTFILE="${PROJECT_BINARY_DIR}/autodoc/doxygen_config"
            -DPROJECT_NAME=${PROJECT_NAME}
            -Dpackage=${PROJECT_NAME}
            -DPROJECT_NUMBER="${DRACO_VERSION_FULL}"
            -Dnumber="${DRACO_VERSION_FULL}"
            -DOUTPUT_DIRECTORY="${DOXYGEN_OUTPUT_DIR}"
            -Ddoxygen_output_top="${DOXYGEN_OUTPUT_DIR}"
            -DINPUT=${DOXYGEN_INPUT}
            -Ddoxygen_input=${DOXYGEN_INPUT}
            -DEXAMPLE_PATH=${DOXYGEN_EXAMPLES}
            -Ddoxygen_examples=${DOXYGEN_EXAMPLES}
            -DHTML_OUTPUT=${DOXYGEN_HTML_OUTPUT}
            -Ddoxygen_html_output=${DOXYGEN_HTML_OUTPUT}
            -DDOTFILE_DIRS=${DOXYGEN_DOTFILE_DIRS}
            -Ddotpath=${DOXYGEN_DOTFILE_DIRS}
            -Dlatex_yes_no=NO
            -Ddoxygen_latex_output=${DOXYGEN_HTML_OUTPUT}/../latex
            -Dhtmlheader=${htmlheader}
            -Dhtmlfooter=${htmlfooter}
            -P "${PROJECT_SOURCE_DIR}/config/configureFileOnMake.cmake"
    DEPENDS "${PROJECT_SOURCE_DIR}/config/doxygen_config.in"
    )
  configure_file(
    ${PROJECT_SOURCE_DIR}/autodoc/html/footer.html.in
    ${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/footer.html
    @ONLY
    )
  configure_file(
    ${PROJECT_SOURCE_DIR}/autodoc/html/header.html.in
    ${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/header.html
    @ONLY
    )
  configure_file(
    ${PROJECT_SOURCE_DIR}/autodoc/html/draco.css
    ${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/draco.css
    COPYONLY
    )
  add_custom_command(
    OUTPUT  "${PROJECT_BINARY_DIR}/autodoc/mainpage.dcc"
    COMMAND "${CMAKE_COMMAND}" 
            -DINFILE="${PROJECT_SOURCE_DIR}/autodoc/mainpage.dcc.in"
            -DOUTFILE="${PROJECT_BINARY_DIR}/autodoc/mainpage.dcc"
            -DCOMP_LINKS=${COMP_LINKS}
            -DPACKAGE_LINKS=${PACKAGE_LINKS}
            -P "${PROJECT_SOURCE_DIR}/config/configureFileOnMake.cmake"
    DEPENDS "${PROJECT_SOURCE_DIR}/autodoc/mainpage.dcc.in"
    COMMENT "Doxygen docs will be installed at"
            "${DOXYGEN_HTML_OUTPUT}.  You must run 'make install'"
  )

  # Copy image files
  # 
  # list of source files
  file( GLOB DOXYGEN_IMAGE_SOURCE_FILES *.jpg *.png html/*.jpg)
  # placeholder for final copy of images.
  set( DOXYGEN_IMAGE_FILES "")
  # build a list of images with full path
  # copy the files from source to binary directory.
  foreach( file ${DOXYGEN_IMAGE_SOURCE_FILES} )
    get_filename_component(tmp ${file} NAME)
    set( DOXYGEN_IMAGE_FILES 
      "${DOXYGEN_IMAGE_FILES} ${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/${tmp}")
    configure_file(
      "${file}"
      "${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/${tmp}"
      COPYONLY
      )
  endforeach( file ${DOXYGEN_IMAGE_SOURCE_FILES} )

  # The actual doxygen command that will be run
  add_custom_command(
    OUTPUT 
       ${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/index.html
       ${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/always_out_of_date.html
    COMMAND "${DOXYGEN_EXECUTABLE}"
    ARGS    doxygen_config
    DEPENDS
      ${DOXYGEN_OUTPUT_DIR}/doxygen_config
      ${DOXYGEN_OUTPUT_DIR}/mainpage.dcc
    WORKING_DIRECTORY ${DOXYGEN_OUTPUT_DIR}
    )

  # Target for building doxygen files.
  add_custom_target( autodoc
    # This DEPENDS links the 'autodoc' make target to the custom command above.    
    # DEPENDS ${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/index.html
    DEPENDS ${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/always_out_of_date.html
    COMMENT "Building Doxygen based autodocs (HTML)...  After running 'make install' your docs will be at ${CMAKE_INSTALL_PREFIX}/autodoc/${DOXYGEN_HTML_OUTPUT}/index.html."
    )

  install( 
    DIRECTORY ${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}
    DESTINATION autodoc
    )

  install(
    FILES ${DOXYGEN_OUTPUT_DIR}/Draco.tag
    DESTINATION autodoc
    OPTIONAL
    )

endif( DOXYGEN_FOUND )
  
# To Do:
# 1.  How to 'make clean'?

