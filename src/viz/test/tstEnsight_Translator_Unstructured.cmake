#-----------------------------*-cmake-*----------------------------------------#
# file   viz/test/tstEnsight_Translator_Unstructured.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   Friday, Aug 28, 2015, 14:27 pm
# brief  This is a CTest script that is used to check the output from
#        viz/test/tstEnsight_Translator_Unstructured
# note   Copyright (C) 2016, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id: CMakeLists.txt 6721 2012-08-30 20:38:59Z gaber $
#------------------------------------------------------------------------------#

# See examples at config/ApplicationUnitTest.cmake and diagnostics/test/tDracoInfo.cmake

#------------------------------------------------------------------------------#
# Setup the CMake based ApplicationUnitTest environment
set( CMAKE_MODULE_PATH ${DRACO_CONFIG_DIR} )
include( ApplicationUnitTest )

# Setup and Sanity check provides:
aut_setup()

##---------------------------------------------------------------------------##
# Run the application and capture the output.
# Variables available for inspection:
#   ${testres} contains the return code
#   ${testout} contains stdout
#   ${testerror} contains stderr
aut_runTests()

##---------------------------------------------------------------------------##
## ASCII files: Run numdiff on the following files:
##---------------------------------------------------------------------------##
set( vars
  geo
  Temperatures
  Pressure
  Velocity
  Densities
  )

# testproblem_ensight directory
foreach( var ${vars} )
  aut_numdiff_2files(
    ${PROJECT_BINARY_DIR}/unstr2d_testproblem_ensight/${var}/data.0001
    ${PROJECT_SOURCE_DIR}/unstr2d_bench/${var}.0001 )
endforeach()

##---------------------------------------------------------------------------##
## Binary files: Run numdiff on the following files:
## - only for Linux.
##---------------------------------------------------------------------------##
macro( aut_diff_2files file1 file2 )

  if( NOT EXISTS ${file1} )
    message( FATAL_ERROR "Specified file1 = ${file1} does not exist." )
  endif()
  if( NOT EXISTS ${file2} )
    message( FATAL_ERROR "Specified file2 = ${file2} does not exist." )
  endif()

  # Assume additional arguments are to be passed to diff

  find_program( exediff diff )
  if( NOT EXISTS ${exediff} )
    message( FATAL_ERROR "Diff not found in PATH")
  endif()
  if( VERBOSE_DEBUG )
    message("   exediff = ${exediff}" )
  endif()

  message("Comparing files:
${exediff} ${ARGV2} ${ARGV3} ${ARGV4} ${ARGV5} ${ARGV6}
   ${file1} \\
   ${file2}")
  execute_process(
    COMMAND ${exediff} ${ARGV2} ${ARGV3} ${ARGV4} ${ARGV5} ${ARGV6} ${file1} ${file2}
    RESULT_VARIABLE diffres
    OUTPUT_VARIABLE diffout
    ERROR_VARIABLE  differror
    )
  if( ${diffres} STREQUAL 0 )
    PASSMSG("files match!
")
  else()
    FAILMSG("files do not match.
diff output = ${diffout}" )
  endif()
endmacro()

# Comparing binary files is expected to fail on powerpc architectures
# because the Endianess of the machine is different than for the
# machine where the gold file was created.  So, for powerpc, skip
# these comparisons.

set( little-endian TRUE )
if( "${CMAKE_HOST_SYSTEM_PROCESSOR}" STREQUAL "powerpc64" OR
    "${CMAKE_HOST_SYSTEM_PROCESSOR}" STREQUAL "ppc64" )
  set( little-endian FALSE )
endif()

if( little-endian )

  # testproblem_binary_ensight directory
  foreach( var ${vars} )
    aut_diff_2files(
      ${PROJECT_BINARY_DIR}/unstr2d_testproblem_binary_ensight/${var}/data.0001
      ${PROJECT_SOURCE_DIR}/unstr2d_bench/${var}.bin.0001 )
  endforeach()

endif()

##---------------------------------------------------------------------------##
## Final report
##---------------------------------------------------------------------------##
aut_report()

##---------------------------------------------------------------------------##
## End
##---------------------------------------------------------------------------##
