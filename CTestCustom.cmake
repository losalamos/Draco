#-----------------------------*-cmake-*----------------------------------------#
# file   CTestCustom.cmake
# brief  Custom configuration for CTest/CDash.
# note   Copyright (C) 2010-2012 LANS, LLC
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

# This file must be found in the root build tree.
# http://www.vtk.org/Wiki/CMake_Testing_With_CTest
# https://www.rad.upenn.edu/sbia/software/doxygen/basis/trunk/html/CTestCustom_8cmake_source.html


##---------------------------------------------------------------------------##
## Warnings
##---------------------------------------------------------------------------##

# Extra matches for warnings:
# set( CTEST_CUSTOM_WARNING_MATCH
#   ${CTEST_CUSTOM_WARNING_MATCH}
#   "{standard input}:[0-9][0-9]*: [Ww]arning: "
#   "{standard input}:[0-9][0-9]*: WARNING: "
#   )

# specialization for machines
# if( "@CMAKE_SYSTEM@" MATCHES "OSF" )
#   set( CTEST_CUSTOM_WARNING_EXCEPTION
#     ${CTEST_CUSTOM_WARNING_EXCEPTION}
#     "XdmfDOM"
#     "XdmfExpr"
#     )
# endif()

# Exceptions 
# set( CTEST_CUSTOM_WARNING_EXCEPTION
#   ${CTEST_CUSTOM_WARNING_EXCEPTION}
#   "tcl8.4.5/[^/]+/../[^/]+/[^.]+.c[:\"]"
#   "Utilities/vtkmpeg2/"
#   "warning LNK44221blahblah"
#   "myspecial/path/to/something/"
#   "myvendorexception"
#   )

# specialization for machines
if( APPLE )
   set( CTEST_CUSTOM_WARNING_EXCEPTION
      ${CTEST_CUSTOM_WARNING_EXCEPTION}
      "has no symbols"
      )
endif()


# specify maximum number of warnings to display
#set( CTEST_CUSTOM_MAXIMUM_NUMBER_OF_WARNINGS "100" )

##---------------------------------------------------------------------------##
## Errors
##---------------------------------------------------------------------------##

## @brief Match expressions for error messages.
# set( CTEST_CUSTOM_ERROR_MATCH
#    ${CTEST_CUSTOM_ERROR_MATCH} # keep current error matches
#    "[0-9][0-9]*: ERROR "       # add match expressions on separate lines
#    "[0-9][0-9]*: [Ee]rror "
#    )

## @brief Match expressions for ignored error messages.
# set( CTEST_CUSTOM_ERROR_EXCEPTION
#    ${CTEST_CUSTOM_ERROR_EXCEPTION} # keep current error exceptions
#    #   "ExampleExec-1.0"           # add exception expressions on separate lines
#    )

# specify maximum number of errors to display
#set( CTEST_CUSTOM_MAXIMUM_NUMBER_OF_ERRORS "100" )

##---------------------------------------------------------------------------##
## Tests
##---------------------------------------------------------------------------##

# ## @brief Specify tests which should be ignored during the test stage.
# set( CTEST_CUSTOM_TESTS_IGNORE ${CTEST_CUSTOM_TESTS_IGNORE} "" )

# ## @brief Specify command to execute before execution of any test during test stage.
# set( CTEST_CUSTOM_PRE_TEST ${CTEST_CUSTOM_PRE_TEST} "" )

# ## @brief Specify command to execute at the end of the test stage.
# set( CTEST_CUSTOM_POST_TEST ${CTEST_CUSTOM_POST_TEST} "" )



##---------------------------------------------------------------------------##
# Code Coverage and Dynamic Analysis Settings
#
# http://www.vtk.org/Wiki/CMake_Testing_With_CTest#Dynamic_Analysis
##---------------------------------------------------------------------------##


# What tool should we use (this should be in your CMakeCache.txt):
# MEMORYCHECK_COMMAND:FILEPATH=/home/kitware/local/bin/valgrind
# PURIFYCOMMAND:FILEPATH=c:/Progra~1/Rational/common/purify.exe

# Add extra options by specifying MEMORYCHECK_COMMAND_OPTIONS and
# MEMORYCHECK_SUPPRESSIONS_FILE.  
set( MEMORYCHECK_SUPPRESSIONS_FILE
   "${CTEST_SCRIPT_DIRECTORY}/valgrind_suppress.txt"
   CACHE FILEPATH
  "File that contains suppressions for the memory checker" )

# Files for exclusion:
set( CTEST_CUSTOM_MEMCHECK_IGNORE
  ${CTEST_CUSTOM_MEMCHECK_IGNORE}
#     test1
#     tstbubbagump 
)

# CTEST_CUSTOM_COVERAGE_EXCLUDE is a list of regular expressions. Any
# file name that matches any of the regular expressions in the list is
# excluded from the reported coverage data.
set( CTEST_CUSTOM_COVERAGE_EXCLUDE 
  ${CTEST_CUSTOM_COVERAGE_EXCLUDE}

  # don't report on actual unit tests
  "/tests/"
  "tests.cpp"
  "tests/tst*.cpp"
  "/src/pkg/tests/tstXercesConfig.cpp"
  )

## @brief Specify additional files which should be considered for coverage report.
#
# Note that the expressions here are globbing expression as
# interpreted by CMake's file(GLOB) command, not regular expressions.
set( CTEST_EXTRA_COVERAGE_GLOB ${CTEST_EXTRA_COVERAGE_GLOB} )
foreach( extension IN ITEMS cc hh )
  list (APPEND CTEST_EXTRA_COVERAGE_GLOB
    "${PROJECT_SOURCE_DIR}/src/*\\.${extension}" )
endforeach()
