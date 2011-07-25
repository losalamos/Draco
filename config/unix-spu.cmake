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
# Compiler Flags
# 

# Special inlining is required on the SPU.  See
# http://gcc.gnu.org/onlinedocs/gcc-4.1.1/gcc/Optimize-Options.html
# for discussion
#
# -finline-limit=n will set:
#    max-inline-insns-single -> n/2
#    max-inline-insns-auto   -> n/2
#    min-inline-insns        -> min(130,n/4)
#    max-inline-insns-rtl    -> n

if( CMAKE_GENERATOR STREQUAL "Unix Makefiles" )

   set( CMAKE_C_FLAGS                "-fstack-check -W -Wall -Winline -fno-exceptions -fno-rtti -march=celledp -DADDRESSING_64 -D__USING_GCC -DEDP -DCACHE_LINE_SIZE=128 -include spu_intrinsics.h" )

   set( CMAKE_C_FLAGS_DEBUG          "-O0 -gdwarf-2 -finline-limit=100 --param large-function-growth=100" )
   set( CMAKE_C_FLAGS_RELEASE        "-O3 -finline-limit=1500 --param large-function-growth=2850" )
   # 2011/07/25 (kt): Looks like we need max-inline-insns-single>800, large-function-growth>2750
   # Setting max-inline-insns-single to a larger might eliminate
   # additional warnings about failed inlining, however doing so causes
   # gcc-4.1.1 to crash
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
