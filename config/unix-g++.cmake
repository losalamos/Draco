#-----------------------------*-cmake-*----------------------------------------#
# file   config/unix-g++.cmake
# author Kelly Thompson 
# date   2008 May 30
# brief  Establish flags for Windows - MSVC
# note   Copyright © 2010 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

#
# Sanity Checks
# 

if( NOT CMAKE_COMPILER_IS_GNUCC )
  message( FATAL_ERROR "If CC is not GNUCC, then we shouldn't have ended up here.  Something is really wrong with the build system. " )
endif( NOT CMAKE_COMPILER_IS_GNUCC )

#
# C++ libraries required by Fortran linker
# 

execute_process( 
  COMMAND ${CMAKE_C_COMPILER} -print-libgcc-file-name
  TIMEOUT 5
  RESULT_VARIABLE tmp
  OUTPUT_VARIABLE libgcc_path
  ERROR_VARIABLE err
  )
get_filename_component( libgcc_path ${libgcc_path} PATH )
execute_process( 
  COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=libstdc++.so
  TIMEOUT 5
  RESULT_VARIABLE tmp
  OUTPUT_VARIABLE libstdcpp_so_loc
  ERROR_VARIABLE err
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )
get_filename_component( libstdcpp_so_loc ${libstdcpp_so_loc} ABSOLUTE )
execute_process( 
  COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=libgcc_s.so
  TIMEOUT 5
  RESULT_VARIABLE tmp
  OUTPUT_VARIABLE libgcc_s_so_loc
  ERROR_VARIABLE err
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )
get_filename_component( libgcc_s_so_loc ${libgcc_s_so_loc} ABSOLUTE )
set( GCC_LIBRARIES 
  ${libstdcpp_so_loc}
  ${libgcc_s_so_loc}
  )
#message(   "   - GNU C++  : ${libstdcpp_so_loc}" )
#message(   "   -          : ${libgcc_s_so_loc}" )

#
# config.h settings
#

execute_process(
  COMMAND ${CMAKE_C_COMPILER} --version
  OUTPUT_VARIABLE ABS_C_COMPILER_VER
  )
string( REGEX REPLACE "Copyright.*" " " 
  ABS_C_COMPILER_VER ${ABS_C_COMPILER_VER} )
string( STRIP ${ABS_C_COMPILER_VER} ABS_C_COMPILER_VER )

execute_process(
  COMMAND ${CMAKE_CXX_COMPILER} --version
  OUTPUT_VARIABLE ABS_CXX_COMPILER_VER
  )
string( REGEX REPLACE "Copyright.*" " " 
  ABS_CXX_COMPILER_VER ${ABS_CXX_COMPILER_VER} )
string( STRIP ${ABS_CXX_COMPILER_VER} ABS_CXX_COMPILER_VER )


#
# Compiler Flags
# 

# Flags from Draco autoconf build system:
# -ansi -pedantic
# -Wnon-virtual-dtor 
# -Wreturn-type 
# -Wno-long-long
# -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC
# -Wextra
# -Weffc++

IF( CMAKE_GENERATOR STREQUAL "Unix Makefiles" )
  set( DRACO_C_FLAGS                "-fPIC -Wcast-align -Wpointer-arith -Wall" )
  set( DRACO_C_FLAGS_DEBUG          "-g -fno-inline -fno-eliminate-unused-debug-types -O0 -Wextra -DDEBUG")
  set( DRACO_C_FLAGS_RELEASE        "-O3 -funroll-loops -march=k8 -DNDEBUG" )
  set( DRACO_C_FLAGS_MINSIZEREL     "${DRACO_C_FLAGS_RELEASE}" )
  set( DRACO_C_FLAGS_RELWITHDEBINFO " -g -fno-inline -fno-eliminate-unused-debug-types -O0 -Wextra -O3 -funroll-loops -march=k8" )

  set( DRACO_CXX_FLAGS                "${DRACO_C_FLAGS}" )
  set( DRACO_CXX_FLAGS_DEBUG          "${DRACO_C_FLAGS_DEBUG} -ansi -pedantic -Woverloaded-virtual -Wno-long-long")
  set( DRACO_CXX_FLAGS_RELEASE        "${DRACO_C_FLAGS_RELEASE}")
  set( DRACO_CXX_FLAGS_MINSIZEREL     "${DRACO_CXX_FLAGS_RELEASE}")
  set( DRACO_CXX_FLAGS_RELWITHDEBINFO "${DRACO_C_FLAGS_RELWITHDEBINFO}" )
ENDIF( CMAKE_GENERATOR STREQUAL "Unix Makefiles" )


string( TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UPPER )
if( ${CMAKE_BUILD_TYPE_UPPER} MATCHES "DEBUG" )
   option( GCC_ENABLE_ALL_WARNINGS 
      "Add \"-Weffc++\" to the compile options (only available for DEBUG builds)." OFF )
   option( GCC_ENABLE_GLIBCXX_DEBUG "Use special version of libc.so that includes STL bounds checking (only available for DEBUG builds)." OFF )
   if( GCC_ENABLE_ALL_WARNINGS )
      set( DRACO_CXX_FLAGS_DEBUG "${DRACO_CXX_FLAGS_DEBUG} -Weffc++" )
      # Force update the CMAKE_CXX_FLAGS (see bottom of this file)
      set( CXX_FLAGS_INITIALIZED "" CACHE INTERNAL "using draco settings." FORCE )
   endif()
   if( GCC_ENABLE_GLIBCXX_DEBUG )
      set( DRACO_CXX_FLAGS_DEBUG "${DRACO_CXX_FLAGS_DEBUG} -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC" )
      # Force update the CMAKE_CXX_FLAGS (see bottom of this file)
      set( CXX_FLAGS_INITIALIZED "" CACHE INTERNAL "using draco settings." FORCE )
   endif()
endif()

if( ENABLE_OPENMP )
  set( DRACO_C_FLAGS   "${DRACO_C_FLAGS} -fopenmp" )
  set( DRACO_CXX_FLAGS "${DRACO_CXX_FLAGS} -fopenmp" )

  # When compiling F90 that links in C++-based libraries, we will need
  # librt added to the link line.
  execute_process( 
    COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=librt.so
    TIMEOUT 5
    RESULT_VARIABLE tmp
    OUTPUT_VARIABLE librt_so_loc
    ERROR_VARIABLE err
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  get_filename_component( librt_so_loc ${librt_so_loc} ABSOLUTE )
  set( GCC_LIBRARIES ${GCC_LIBRARIES} ${librt_so_loc} )
endif()

option( ENABLE_C_CODECOVERAGE "Instrument for C/C++ code coverage analysis?" OFF )
if( ENABLE_C_CODECOVERAGE )
  find_program( COVERAGE_COMMAND gcov )
  set( DRACO_C_FLAGS_DEBUG     "${DRACO_C_FLAGS_DEBUG} -O0 -fprofile-arcs -ftest-coverage" )
  set( DRACO_CXX_FLAGS_DEBUG   "${DRACO_C_FLAGS_DEBUG}")
  set( CMAKE_LDFLAGS           "-fprofile-arcs -ftest-coverage" )

  # When compiling F90 that links in C++-based libraries, we will need
  # libgcov added to the link line.
  execute_process( 
    COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=libgcov.a
    TIMEOUT 5
    RESULT_VARIABLE tmp
    OUTPUT_VARIABLE libgcov_a_loc
    ERROR_VARIABLE err
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  get_filename_component( libgcov_a_loc ${libgcov_a_loc} ABSOLUTE )
  set( GCC_LIBRARIES ${GCC_LIBRARIES} ${libgcov_a_loc} )
endif( ENABLE_C_CODECOVERAGE )

# Save the Draco default values to the cache file
if( "${CXX_FLAGS_INITIALIZED}no" STREQUAL "no" )
   set( CXX_FLAGS_INITIALIZED "yes" CACHE INTERNAL "using draco settings." )
   set( CMAKE_C_FLAGS                "${DRACO_C_FLAGS}"                CACHE STRING "compiler flags" FORCE )
   set( CMAKE_C_FLAGS_DEBUG          "${DRACO_C_FLAGS_DEBUG}"          CACHE STRING "compiler flags" FORCE ) 
   set( CMAKE_C_FLAGS_RELEASE        "${DRACO_C_FLAGS_RELEASE}"        CACHE STRING "compiler flags" FORCE )
   set( CMAKE_C_FLAGS_MINSIZEREL     "${DRACO_C_FLAGS_MINSIZEREL}"     CACHE STRING "compiler flags" FORCE )
   set( CMAKE_C_FLAGS_RELWITHDEBINFO "${DRACO_C_FLAGS_RELWITHDEBINFO}" CACHE STRING "compiler flags" FORCE )
   set( CMAKE_CXX_FLAGS                "${DRACO_CXX_FLAGS}"                CACHE STRING "compiler flags" FORCE )
   set( CMAKE_CXX_FLAGS_DEBUG          "${DRACO_CXX_FLAGS_DEBUG}"          CACHE STRING "compiler flags" FORCE ) 
   set( CMAKE_CXX_FLAGS_RELEASE        "${DRACO_CXX_FLAGS_RELEASE}"        CACHE STRING "compiler flags" FORCE )
   set( CMAKE_CXX_FLAGS_MINSIZEREL     "${DRACO_CXX_FLAGS_MINSIZEREL}"     CACHE STRING "compiler flags" FORCE )
   set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${DRACO_CXX_FLAGS_RELWITHDEBINFO}" CACHE STRING "compiler flags" FORCE )
endif()



#------------------------------------------------------------------------------#
# End config/unix-g++.cmake
#------------------------------------------------------------------------------#
