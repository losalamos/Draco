#-----------------------------*-cmake-*----------------------------------------#
# file   config/ApplicationUnitTest.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   Monday, Nov 19, 2012, 16:21 pm
# brief  Provide macros that aid in creating unit tests that run
#        interactive user codes (i.e.: run a binary that reads an
#        input file and diff the resulting output file).
# note   Copyright (C) 2012, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id: CMakeLists.txt 6732 2012-09-05 22:28:18Z kellyt $
#------------------------------------------------------------------------------#

##---------------------------------------------------------------------------##
## Check values for $APP
##
## Returns:
##   APP       - cleaned up path to executable
##   STDINFILE - cleaned up path to intput file.
##   BINDIR    - Directory location of binary file
##   PROJECT_BINARY_DIR - Parent directory of BINDIR
##   GOLDFILE  - cleaned up path to gold standard file.
##   OUTFILE   - Output filename derived from GOLDFILE.
##
## if VERBOSE is set, also echo input values.
##   APP       - path name for executable.
##   STDINFILE - optional input file
##   GOLDFILE  - optional gold standard file.
##---------------------------------------------------------------------------##
macro( aut_setup)
   
   if( VERBOSE )
      message("Running tQueryEospac.cmake with the following parameters:")
      message("   APP       = ${APP}")
      if( STDINFILE )
         message("   STDINFILE = ${STDINFILE}")
      endif()
      if( GOLDFILE )
         message("   GOLDFILE = ${GOLDFILE}" )
      endif()
   endif()

   # Setup and sanity check 

   if( "${APP}x" STREQUAL "x" )
      message( FATAL_ERROR "You must provide a value for APP." )
   endif()
   
   get_filename_component( APP ${APP} ABSOLUTE )
   if( STDINFILE )
      get_filename_component( STDINFILE ${STDINFILE} ABSOLUTE )
   endif()
   get_filename_component( BINDIR ${APP} PATH )
   get_filename_component( PROJECT_BINARY_DIR ${BINDIR} PATH )
   if( GOLDFILE )
      get_filename_component( OUTFILE ${GOLDFILE} NAME_WE )
   else()
      get_filename_component( OUTFILE ${APP} NAME_WE )
   endif()
   set( OUTFILE "${CMAKE_CURRENT_BINARY_DIR}/${OUTFILE}.out")
   
   message("Testing ${APP}")
   
   if( NOT EXISTS "${APP}" )
      message( FATAL_ERROR "Cannot find ${APP}")
   endif()

   set( numpasses 0 )
   set( numfails  0 )

endmacro()

##---------------------------------------------------------------------------##
## PASSMSG/FAILMSG
##---------------------------------------------------------------------------##

function(PASSMSG msg)
    math( EXPR numpasses "${numpasses} + 1" )
    message( "Test Passes: ${msg}")
endfunction()

function(FAILMSG msg)
    math( EXPR numfails "${numfails} + 1" )
    message( "Test Fails: ${msg}")
endfunction()
