#-----------------------------*-cmake-*----------------------------------------#
# file   config/unix-spu.cmake
# author Kelly Thompson 
# date   2011 May 11
# brief  Establish flags for Roadrunner (PowerPC)
# note   Copyright © 2011 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

#
# Sanity Checks
# 

#
# C++ libraries required by Fortran linker
# 

# execute_process( 
#   COMMAND ${CMAKE_C_COMPILER} -print-libgcc-file-name
#   TIMEOUT 5
#   RESULT_VARIABLE tmp
#   OUTPUT_VARIABLE libgcc_path
#   ERROR_VARIABLE err
#   )
# get_filename_component( libgcc_path ${libgcc_path} PATH )
# execute_process( 
#   COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=libstdc++.so
#   TIMEOUT 5
#   RESULT_VARIABLE tmp
#   OUTPUT_VARIABLE libstdcpp_so_loc
#   ERROR_VARIABLE err
#   OUTPUT_STRIP_TRAILING_WHITESPACE
#   )
# get_filename_component( libstdcpp_so_loc ${libstdcpp_so_loc} ABSOLUTE )
# execute_process( 
#   COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=libgcc_s.so
#   TIMEOUT 5
#   RESULT_VARIABLE tmp
#   OUTPUT_VARIABLE libgcc_s_so_loc
#   ERROR_VARIABLE err
#   OUTPUT_STRIP_TRAILING_WHITESPACE
#   )
# get_filename_component( libgcc_s_so_loc ${libgcc_s_so_loc} ABSOLUTE )
# set( GCC_LIBRARIES 
#   ${libstdcpp_so_loc}
#   ${libgcc_s_so_loc}
#   )
#message(   "   - GNU C++  : ${libstdcpp_so_loc}" )
#message(   "   -          : ${libgcc_s_so_loc}" )

#
# config.h settings
#

# execute_process(
#   COMMAND ${CMAKE_C_COMPILER} --version
#   OUTPUT_VARIABLE ABS_C_COMPILER_VER
#   )
# string( REGEX REPLACE "Copyright.*" " " 
#   ABS_C_COMPILER_VER ${ABS_C_COMPILER_VER} )
# string( STRIP ${ABS_C_COMPILER_VER} ABS_C_COMPILER_VER )

# execute_process(
#   COMMAND ${CMAKE_CXX_COMPILER} --version
#   OUTPUT_VARIABLE ABS_CXX_COMPILER_VER
#   )
# string( REGEX REPLACE "Copyright.*" " " 
#   ABS_CXX_COMPILER_VER ${ABS_CXX_COMPILER_VER} )
# string( STRIP ${ABS_CXX_COMPILER_VER} ABS_CXX_COMPILER_VER )


#
# Compiler Flags
# 

if( CMAKE_GENERATOR STREQUAL "Unix Makefiles" )

   set( CMAKE_C_FLAGS                "-fstack-check -W -Wall -Winline -fno-exceptions -fno-rtti -march=celledp -DADDRESSING_64 -D__USING_GCC -DEDP -DCACHE_LINE_SIZE=128 -include spu_intrinsics.h" )

   set( CMAKE_C_FLAGS_DEBUG          "-O0 -gdwarf-2 -finline-limit=100 --param large-function-growth=100" )
   set( CMAKE_C_FLAGS_RELEASE        "-O3 -finline-limit=1500 --param large-function-growth=1700" )
   set( CMAKE_C_FLAGS_MINSIZEREL     "-O3" )
   set( CMAKE_C_FLAGS_RELWITHDEBINFO "-O3 -gdwarf-2 -finline-limit=100 --param large-function-growth=100" )
   
   set( CMAKE_CXX_FLAGS                "${CMAKE_C_FLAGS}" )
   set( CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_C_FLAGS_DEBUG}")
   set( CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_C_FLAGS_RELEASE}")
   set( CMAKE_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_RELEASE}")
   set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}"
      )
   
   # remove -rdynamic from link line when creating an executable
   set( CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS )

   SET(CMAKE_C_ARCHIVE_CREATE "<CMAKE_AR> <LINK_FLAGS> -qcs <TARGET> <OBJECTS>")
   #SET(CMAKE_C_ARCHIVE_APPEND "<CMAKE_AR> <LINK_FLAGS> r <TARGET> <OBJECTS>")
   SET(CMAKE_CXX_ARCHIVE_CREATE "<CMAKE_AR> <LINK_FLAGS> -qcs <TARGET> <OBJECTS>")
   #SET(CMAKE_CXX_ARCHIVE_APPEND "<CMAKE_AR> <LINK_FLAGS> r <TARGET>
   #<OBJECTS>")

else()
   message( FATAL_ERROR "
I dont' know how to setup the spu-g++ compiler for build systems other 
than Unix Makefiles.")
endif()

#------------------------------------------------------------------------------#
# End config/unix-g++.cmake
#------------------------------------------------------------------------------#
