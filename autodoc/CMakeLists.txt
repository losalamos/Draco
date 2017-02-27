#-----------------------------*-cmake-*----------------------------------------#
# file   draco/autodoc/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 Oct 14
# brief  Instructions for building Doxygen documentation
# note   Copyright (C) 2010-2017 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

# Doxygen
find_package( Doxygen QUIET )

if( DOXYGEN_FOUND )

   # if AUTODOCDIR is set in environment (or make command line), create a CMake
   # variable with this value.
   if( NOT "$ENV{AUTODOCDIR}x" STREQUAL "x" )
      set( AUTODOCDIR "$ENV{AUTODOCDIR}" )
   endif()
   # if AUTODOCDIR is set, use it, otherwise, default to CMAKE_INSTALL_PREFIX.
   if( "${AUTODOCDIR}x" STREQUAL "x" )
      set( DOXYGEN_OUTPUT_DIR ${CMAKE_INSTALL_PREFIX}/autodoc )
   else()
      set( DOXYGEN_OUTPUT_DIR "${AUTODOCDIR}" )
   endif()

  # Source directories that should be examined by doxygen Sets
  # ${DOXYGEN_INPUT}.  This value is embedded in the doxygen_config
  # file below.
  set( DOXYGEN_INPUT
     "${PROJECT_SOURCE_DIR}/autodoc ${PROJECT_BINARY_DIR}/autodoc" )

  file( GLOB package_list ${PROJECT_SOURCE_DIR}/src/* )
  foreach( package ${package_list} )
     if( EXISTS ${package}/CMakeLists.txt )
        set( DOXYGEN_INPUT "${DOXYGEN_INPUT} ${package}" )
        # message( "Doxygen will process ${package}" )
     endif()
     if( EXISTS ${package}/test/CMakeLists.txt )
        set( DOXYGEN_INPUT "${DOXYGEN_INPUT} ${package}/test" )
        set( DOXYGEN_EXAMPLES "${DOXYGEN_EXAMPLES} ${package}/test" )
        # message( "Doxygen will process examples found at ${package}/test" )
     endif()
     if( EXISTS ${package}/autodoc )
        set( DOXYGEN_INPUT "${DOXYGEN_INPUT} ${package}/autodoc" )
        # message( "Doxygen will process ${package}/autodoc" )
     endif()
  endforeach()

  # Also look for files that have been configured (.in files) and
  # placed in the BINARY_DIR.
  file( GLOB package_list ${PROJECT_BINARY_DIR}/src/* )
  foreach( package ${package_list} )
     # message( "Looking at ${package}")
     # pick up processed .dcc files
     if( EXISTS ${package}/autodoc )
        set( DOXYGEN_INPUT "${DOXYGEN_INPUT} ${package}/autodoc" )
        # message( "Doxygen will process ${package}/autodoc" )
     endif()
  endforeach()

  set( DOXYGEN_HTML_OUTPUT draco )
  set( DOXYGEN_DOTFILE_DIRS dot )

  set( htmlfooter "${DOXYGEN_HTML_OUTPUT}/footer.html")
  set( htmlheader "${DOXYGEN_HTML_OUTPUT}/header.html")
  set( project_brief "Object-oriented components for numerically intensive, parallel, radiation transport applications.")
  string( REPLACE " " "___" project_brief ${project_brief} )

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
            -Dproject_brief=\"${project_brief}\"
            -Dlatex_yes_no=NO
            -Ddoxygen_latex_output=${DOXYGEN_HTML_OUTPUT}/../latex
            -Dhtmlheader=${htmlheader}
            -Dhtmlfooter=${htmlfooter}
            -Dproject_logo=dragon.jpg
            -P "${PROJECT_SOURCE_DIR}/config/configureFileOnMake.cmake"
    DEPENDS ${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/always_out_of_date.html
    # DEPENDS "${PROJECT_SOURCE_DIR}/config/doxygen_config.in"
    )
  configure_file(
    ${PROJECT_SOURCE_DIR}/autodoc/html/footer.html.in
    ${PROJECT_BINARY_DIR}/autodoc/${DOXYGEN_HTML_OUTPUT}/footer.html
    @ONLY
    )
  configure_file(
    ${PROJECT_SOURCE_DIR}/autodoc/html/header.html.in
    ${PROJECT_BINARY_DIR}/autodoc/${DOXYGEN_HTML_OUTPUT}/header.html
    @ONLY
    )
  configure_file(
    ${PROJECT_SOURCE_DIR}/autodoc/html/draco.css
    ${PROJECT_BINARY_DIR}/autodoc/${DOXYGEN_HTML_OUTPUT}/draco.css
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
    COMMENT "Doxygen docs are installed at"
            "${DOXYGEN_HTML_OUTPUT}.  " # You must run 'make install'"
  )

  # Copy image files
  #
  # list of source files
  file( GLOB DOXYGEN_IMAGE_SOURCE_FILES *.jpg *.png html/*.jpg html/*.png)
  # placeholder for final copy of images.
  set( DOXYGEN_IMAGE_FILES "")
  # build a list of images with full path
  # copy the files from source to binary directory.
  foreach( file ${DOXYGEN_IMAGE_SOURCE_FILES} )
    get_filename_component(filename_nopath ${file} NAME)
    set( DOXYGEN_IMAGE_FILES
      "${DOXYGEN_IMAGE_FILES} ${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/${filename_nopath}")
    configure_file(
      "${file}"
      "${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/${filename_nopath}"
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
      ${PROJECT_BINARY_DIR}/autodoc/doxygen_config
      ${PROJECT_BINARY_DIR}/autodoc/mainpage.dcc
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/autodoc
    )

  # Target for building doxygen files.
  add_custom_target( autodoc
    DEPENDS ${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/always_out_of_date.html
    COMMENT "Building Doxygen based autodocs (HTML)...

*** Autodoc HTML pages are at
    file://${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/index.html

" )

  add_custom_target( installandautodoc
     COMMAND ${CMAKE_MAKE_PROGRAM} install
     COMMAND ${CMAKE_MAKE_PROGRAM} autodoc
     WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  )

endif( DOXYGEN_FOUND )

# To Do:
# 1. How to 'make clean'?
# 2. Consider using tagfiles and building each component as an thread
#    to decrease wall clock time needed to build docs.
