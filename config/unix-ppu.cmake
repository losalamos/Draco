#-----------------------------*-cmake-*----------------------------------------#
# file   config/unix-ppu.cmake
# author Kelly Thompson 
# date   2011 April 15
# brief  Establish flags for Roadrunner (PowerPC)
# note   Copyright © 2011 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

#
# Sanity Checks
# 

# Must use static libraries.
set( DRACO_LIBRARY_TYPE "STATIC" CACHE STRING 
   "Keyword for creating new libraries (STATIC or SHARED)."
   FORCE )

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

# Flags from Draco autoconf build system:
# -ansi -pedantic
# -Wnon-virtual-dtor 
# -Wreturn-type 
# -Wno-long-long
# -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC
# -Wextra
# -Weffc++

# /opt/cell/toolchain/bin/ppu-g++ -o Release.o -c -W -Wall
# -Wno-sign-compare -m64 -mabi=altivec -maltivec -O0 -D__PPU__
# -DADDRESSING_64 -DDBC=7 -gdwarf-2

# -DMESH_EVERY_CYCLE 
# -DHOST_ACCEL_DACS 
# -DPPE_WRITE_BUFFER_DIRECT 
# -DPPE_READ_BUFFER_DIRECT 
# -DACCEL_RECV_NONBLOCKING 
# -DACCEL_SEND_NONBLOCKING 
# -DSHORT_SPE_TALLIES 
# -DCACHE_LINE_SIZE=128 
# -Ids++ Release.cc


if( CMAKE_GENERATOR STREQUAL "Unix Makefiles" )

   set( CMAKE_C_FLAGS                "-W -Wall -Wno-sign-compare -m64 -mabi=altivec -maltivec -D__PPU__ -DADDRESSING_64 -DCACHE_LINE_SIZE=128" )
   set( CMAKE_C_FLAGS_DEBUG          "-O0 -gdwarf-2" )
   set( CMAKE_C_FLAGS_RELEASE        "-O3" )
   set( CMAKE_C_FLAGS_MINSIZEREL     "-O3" )
   set( CMAKE_C_FLAGS_RELWITHDEBINFO "-O3 -gdwarf-2" )
   
   set( CMAKE_CXX_FLAGS                "${CMAKE_C_FLAGS}" )
   set( CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_C_FLAGS_DEBUG}")
   set( CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_C_FLAGS_RELEASE}")
   set( CMAKE_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_RELEASE}")
   set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}"
      )
   
   SET(CMAKE_C_ARCHIVE_CREATE "<CMAKE_AR> <LINK_FLAGS> -qcs <TARGET> <OBJECTS>")
   #SET(CMAKE_C_ARCHIVE_APPEND "<CMAKE_AR> <LINK_FLAGS> r <TARGET> <OBJECTS>")
   SET(CMAKE_CXX_ARCHIVE_CREATE "<CMAKE_AR> <LINK_FLAGS> -qcs <TARGET> <OBJECTS>")
   #SET(CMAKE_CXX_ARCHIVE_APPEND "<CMAKE_AR> <LINK_FLAGS> r <TARGET>
   #<OBJECTS>")

else()
   message( FATAL_ERROR "
I dont' know how to setup the ppu-g++ compiler for build systems other 
than Unix Makefiles.")
endif()

#------------------------------------------------------------------------------#
# End config/unix-g++.cmake
#------------------------------------------------------------------------------#
