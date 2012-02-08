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

#
# config.h settings
#

execute_process(
  COMMAND ${CMAKE_C_COMPILER} --version
  OUTPUT_VARIABLE DBS_C_COMPILER_VER
  )
string( REGEX REPLACE "Copyright.*" " " 
  DBS_C_COMPILER_VER ${DBS_C_COMPILER_VER} )
string( STRIP ${DBS_C_COMPILER_VER} DBS_C_COMPILER_VER )

execute_process(
  COMMAND ${CMAKE_CXX_COMPILER} --version
  OUTPUT_VARIABLE DBS_CXX_COMPILER_VER
  )
string( REGEX REPLACE "Copyright.*" " " 
  DBS_CXX_COMPILER_VER ${DBS_CXX_COMPILER_VER} )
string( STRIP ${DBS_CXX_COMPILER_VER} DBS_CXX_COMPILER_VER )

string( REGEX REPLACE ".*([0-9]).([0-9]).([0-9]).*" "\\1"
   DBS_CXX_COMPILER_VER_MAJOR ${DBS_CXX_COMPILER_VER} )
string( REGEX REPLACE ".*([0-9]).([0-9]).([0-9]).*" "\\2"
   DBS_CXX_COMPILER_VER_MINOR ${DBS_CXX_COMPILER_VER} )

# is this bullseye?
execute_process(
  COMMAND ${CMAKE_CXX_COMPILER} --version
  ERROR_VARIABLE DBS_CXX_IS_BULLSEYE
  OUTPUT_VARIABLE tmp
  )

if( ${DBS_CXX_IS_BULLSEYE} MATCHES BullseyeCoverage )
   set( DBS_CXX_IS_BULLSEYE "ON" )
   set( DBS_CXX_IS_BULLSEYE ${DBS_CXX_IS_BULLSEYE} 
      CACHE BOOL "Are we running Bullseye" )
   mark_as_advanced( DBS_CXX_IS_BULLSEYE )
endif()

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

# -march=core2 | corei7 
# -mtune=core2 | corei7 

if( NOT CXX_FLAGS_INITIALIZED )
   set( CXX_FLAGS_INITIALIZED "yes" CACHE INTERNAL "using draco settings." )

   set( CMAKE_C_FLAGS                "-fPIC -Wcast-align -Wpointer-arith -Wall -fopenmp" )
   set( CMAKE_C_FLAGS_DEBUG          "-g -fno-inline -fno-eliminate-unused-debug-types -O0 -Wextra -DDEBUG")
   set( CMAKE_C_FLAGS_RELEASE        "-O3 -funroll-loops -DNDEBUG" )
   set( CMAKE_C_FLAGS_MINSIZEREL     "${CMAKE_C_FLAGS_RELEASE}" )
   set( CMAKE_C_FLAGS_RELWITHDEBINFO "-O3 -g -fno-eliminate-unused-debug-types -Wextra -funroll-loops" )

   set( CMAKE_CXX_FLAGS                "${CMAKE_C_FLAGS}" )
   set( CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_C_FLAGS_DEBUG} -ansi -pedantic -Woverloaded-virtual -Wno-long-long")
   set( CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_C_FLAGS_RELEASE}")
   set( CMAKE_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_RELEASE}")
   set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}" )

   # XCode needs this extra bit of help for OpenMP support
   set( CMAKE_EXE_LINKER_FLAGS "-fopenmp" )
endif()

# Extra flags for gcc-4.6.2+
# -Wsuggest-attribute=const
if( ${DBS_CXX_COMPILER_VER_MAJOR} GREATER 3 )
   if( ${DBS_CXX_COMPILER_VER_MINOR} GREATER 5 )
      set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wsuggest-attribute=const" )
   endif()
endif()


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

##---------------------------------------------------------------------------##
# Ensure cache values always match current selection
##---------------------------------------------------------------------------##
set( CMAKE_C_FLAGS                "${CMAKE_C_FLAGS}"                CACHE STRING "compiler flags" FORCE )
set( CMAKE_C_FLAGS_DEBUG          "${CMAKE_C_FLAGS_DEBUG}"          CACHE STRING "compiler flags" FORCE ) 
set( CMAKE_C_FLAGS_RELEASE        "${CMAKE_C_FLAGS_RELEASE}"        CACHE STRING "compiler flags" FORCE )
set( CMAKE_C_FLAGS_MINSIZEREL     "${CMAKE_C_FLAGS_MINSIZEREL}"     CACHE STRING "compiler flags" FORCE )
set( CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}" CACHE STRING "compiler flags" FORCE )

set( CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS}"                CACHE STRING "compiler flags" FORCE )
set( CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG}"          CACHE STRING "compiler flags" FORCE ) 
set( CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE}"        CACHE STRING "compiler flags" FORCE )
set( CMAKE_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_MINSIZEREL}"     CACHE STRING "compiler flags" FORCE )
set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}" CACHE STRING "compiler flags" FORCE )

#------------------------------------------------------------------------------#
# End config/unix-g++.cmake
#------------------------------------------------------------------------------#
