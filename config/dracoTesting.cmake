#-----------------------------*-cmake-*----------------------------------------#
# file   config/compilerEnv.cmake
# brief  Default CMake build parameters
# note   Copyright (C) 2016-2019 Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

include_guard(GLOBAL)

include( FeatureSummary )
option( BUILD_TESTING "Should we compile the tests?" ON )
add_feature_info( BUILD_TESTING BUILD_TESTING
   "Turn off to prevent the compilation of unit tests (ctest).")

# how many cores on the local system?
if( UNIX )
  if( EXISTS "/proc/cpuinfo" )
    file( READ "/proc/cpuinfo" cpuinfo )
    # convert one big string into a set of strings, one per line
    string( REGEX REPLACE "\n" ";" cpuinfo ${cpuinfo} )
    set( proc_ids "" )
    foreach( line ${cpuinfo} )
       if( ${line} MATCHES "processor" )
          list( APPEND proc_ids ${line} )
       endif()
    endforeach()
    list( LENGTH proc_ids DRACO_NUM_CORES )
    set( MPIEXEC_MAX_NUMPROCS ${DRACO_NUM_CORES} CACHE STRING
       "Number of cores on the local machine." )
  endif()
endif()

# enable ctest funcitons and run ctest in parallel if we have multiple cores.
if( BUILD_TESTING )
  include(CTest)
  enable_testing()
  # by default do not use parallel build flags (e.g.: -j16)
  cmake_host_system_information( RESULT logical_cores QUERY NUMBER_OF_LOGICAL_CORES )
  set( pbuildtestflags "" )
  if( NOT WIN32 ) # stick with scalar builds for now.
     set( pbuildtestflags "-j${logical_cores}" )
  endif()
  if( ${CMAKE_GENERATOR} MATCHES Ninja )
    add_custom_target( check
      COMMAND "${CMAKE_COMMAND}" --build "${Draco_BINARY_DIR}" -- ${pbuildtestflags}
      COMMAND ${CMAKE_CTEST_COMMAND} ${pbuildtestflags} $$(ARGS) )
  else()
    add_custom_target( check
      COMMAND "${CMAKE_COMMAND}" --build "${Draco_BINARY_DIR}" -- ${pbuildtestflags}
      COMMAND ${CMAKE_CTEST_COMMAND} ${pbuildtestflags} $(ARGS) )
  endif()
endif()

#------------------------------------------------------------------------------#
# End dracoTesting.cmake
#------------------------------------------------------------------------------#
