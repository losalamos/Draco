#-----------------------------*-cmake-*----------------------------------------#
# file   draco/config/configureFileOnMake.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 Oct 14
# brief  Configure a file during the 'make' process instead of during 'cmake'
# note   Copyright (C) 2016-2020 Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

# Use:
#   cmake -PconfigureFileOnMake.cmake -DINFILE=$infile \
#         -DOUTFILE=$outfile -DSUBST_VARIABLE=VALUE ...
#
# Prerequisits:
# - INFILE must be a valid file.
#
# Post
# - OUTFILE will be written (or overwritten).
#
# Suggested use is to generate binary directory files based on
# build-time changes to corresonding source tree .in file. For
# example,
#   add_custom_command(
#    OUTPUT  "${PROJECT_BINARY_DIR}/autodoc/mainpage.dcc"
#    COMMAND "${CMAKE_COMMAND}"
#            -DINFILE="${PROJECT_SOURCE_DIR}/autodoc/mainpage.dcc.in"
#            -DOUTFILE="${PROJECT_BINARY_DIR}/autodoc/mainpage.dcc"
#            -DCOMP_LINKS=${COMP_LINKS}
#            -DPACKAGE_LINKS=${PACKAGE_LINKS}
#            -P "${PROJECT_SOURCE_DIR}/config/configureFileOnMake.cmake"
#    DEPENDS "${PROJECT_SOURCE_DIR}/autodoc/mainpage.dcc.in"
#  )
#
# add_custom_target( mytarget
#   DEPENDS "${PROJECT_BINARY_DIR}/autodoc/mainpage.dcc"
#    COMMENT "Building Doxygen mainpage.dcc (HTML)..."
#    )

if( NOT EXISTS ${INFILE} )
  message( FATAL_ERROR "
INFILE and OUTFILE must be set on command line.  For example,
${CMAKE_COMMAND}
  -P ${PROJECT_SOURCE_DIR}/config/configureFileOnMake.cmake
  -DINFILE=${INFILE}
  -DOUTFILE=${OUTFILE}
  -DPROJECT_NAME=${PROJECT_NAME}
  -DPROJECT_NUMBER=${PROJECT_NUMBER}
  -DOUTPUT_DIRECTORY=${OUTPUT_DIRECTORY}
  -DINPUT=${INPUT}
  -DEXAMPLE_PATH=${EXAMPLE_PATH}
  -DHTML_OUTPUT=${HTML_OUTPUT}
  -DTAGFILES=${TAGFILES}
  -DDOTFILE_DIRS=${DOTFILE_DIRS}
" )
endif()

# Decode "---" as " " for variables passed to this function
string( REPLACE "___" " " project_brief "${project_brief}" )

# Ensure we use native path styles
file( TO_NATIVE_PATH "${PROJECT_SOURCE_DIR}" PROJECT_SOURCE_DIR )
file( TO_NATIVE_PATH "${INFILE}" INFILE )
file( TO_NATIVE_PATH "${OUTFILE}" OUTFILE )
file( TO_NATIVE_PATH "${OUTPUT_DIRECTORY}" OUTPUT_DIRECTORY )
file( TO_NATIVE_PATH "${EXAMPLE_PATH}" EXAMPLE_PATH )
file( TO_NATIVE_PATH "${HTML_OUTPUT}" HTML_OUTPUT )
file( TO_NATIVE_PATH "${DOTFILE_DIRS}" DOTFILE_DIRS )
file( TO_NATIVE_PATH "${DRACO_DIR}" DRACO_DIR )
# INPUT is a list of paths
#message("INPUT = ${INPUT}")
#string( REGEX REPLACE "[ ]" ";" INPUT ${INPUT} )
#message("INPUT = ${INPUT}")
#foreach( item ${INPUT} )
#   file( TO_NATIVE_PATH "${item}" input_path )
#   set( tmp_input "${tmp_input} ${input_path}" )
#   message("
#${item} --> ${tmp_input}
#INPUT = ${input_path}
#")
#endforeach()
#set( INPUT ${tmp_input} )
# TAGFILES is a list of paths
if( TAGFILES )
   string( REGEX REPLACE "([ ])" ";" list ${TAGFILES} )
   foreach( item ${list} )
      file( TO_NATIVE_PATH "${item}" input_path )
      if( ${input_path} MATCHES ".tag" )
         set( tmp_tagfiles "${tmp_tagfiles} ${input_path}" )
      else()
         set( tmp_tagfiles "${tmp_tagfiles}=${input_path}" )
      endif()
   endforeach()
   set( TAGFILES ${tmp_tagfiles} )
endif()
configure_file( ${INFILE} ${OUTFILE} @ONLY )
