#-----------------------------*-cmake-*----------------------------------------#
# file   CTestCustom.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 Dec 8
# brief  Custom configuration for CTest/CDash.
# note   © Copyright 2010 LANS, LLC
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

# This file must be found in the root build tree.
# http://www.vtk.org/Wiki/CMake_Testing_With_CTest

# Extra matches for warnings:
# set( CTEST_CUSTOM_WARNING_MATCH
#   ${CTEST_CUSTOM_WARNING_MATCH}
#   "{standard input}:[0-9][0-9]*: Warning: "
#   )

# specialization for machines
#if( "@CMAKE_SYSTEM@" MATCHES "OSF" )
#  set( CTEST_CUSTOM_WARNING_EXCEPTION
#    ${CTEST_CUSTOM_WARNING_EXCEPTION}
#    "XdmfDOM"
#    "XdmfExpr"
#    )
#endif( "@CMAKE_SYSTEM@" MATCHES "OSF" )

# Exceptions 
set( CTEST_CUSTOM_WARNING_EXCEPTION
  ${CTEST_CUSTOM_WARNING_EXCEPTION}
  "tcl8.4.5/[^/]+/../[^/]+/[^.]+.c[:\"]"
  "Utilities/vtkmpeg2/"
  "warning LNK44221blahblah"
  "myspecial/path/to/something/"
  "myvendorexception"
  )

# ----------------------------------------
# Code Coverage and Dynamic Analysis Settings
#
# http://www.vtk.org/Wiki/CMake_Testing_With_CTest#Dynamic_Analysis
# ----------------------------------------

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
