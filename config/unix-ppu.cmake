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
