#-----------------------------*-cmake-*----------------------------------------#
# file   config/apple-clang.cmake
# brief  Establish flags for Apple OSX
# note   Copyright (C) 2010-2013 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id: unix-g++.cmake 6798 2012-10-09 21:45:34Z kellyt $
#------------------------------------------------------------------------------#

#
# Compiler flag checks
#
include(platform_checks)
query_openmp_availability()

#
# Compiler Flags
#

if( NOT CXX_FLAGS_INITIALIZED )
   set( CXX_FLAGS_INITIALIZED "yes" CACHE INTERNAL "using draco settings." )

   set( CMAKE_C_FLAGS                "-Wcast-align -Wpointer-arith -Wall -Wno-long-long -pedantic" )
   if (NOT ${CMAKE_GENERATOR} MATCHES Xcode AND HAS_MARCH_NATIVE)
      set( CMAKE_C_FLAGS                "${CMAKE_C_FLAGS} -march=native" )
   endif()
   set( CMAKE_C_FLAGS_DEBUG          "-g -fno-inline -O0 -Wextra -DDEBUG")
   set( CMAKE_C_FLAGS_RELEASE        "-O3 -funroll-loops -DNDEBUG" )
   set( CMAKE_C_FLAGS_MINSIZEREL     "${CMAKE_C_FLAGS_RELEASE}" )
   set( CMAKE_C_FLAGS_RELWITHDEBINFO "-O3 -g -Wextra -funroll-loops" )

   set( CMAKE_CXX_FLAGS                "${CMAKE_C_FLAGS}" )
   set( CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_C_FLAGS_DEBUG} -Woverloaded-virtual")
   set( CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_C_FLAGS_RELEASE}")
   set( CMAKE_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_RELEASE}")
   set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}" )

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

# Toggle for C++11 support
toggle_compiler_flag( DRACO_ENABLE_STRICT_ANSI "-std=c++98" "CXX" "")
toggle_compiler_flag( DRACO_ENABLE_STRICT_ANSI "-std=c90"   "C" "")
toggle_compiler_flag( DRACO_ENABLE_C99         "-std=c99" "C" "" )
if( OpenMP_C_FLAGS )
  message("toggle_compiler_flag( OPENMP_FOUND \"${OpenMP_C_FLAGS}\" \"C;CXX;EXE_LINKER\" \"\" )")
  toggle_compiler_flag( OPENMP_FOUND "${OpenMP_C_FLAGS}" "C;CXX;EXE_LINKER" "" )
endif()
# toggle_compiler_flag( DRACO_ENABLE_CXX11 "-stdlib=libc++ -std=c++11" "CXX" "")
toggle_compiler_flag( DRACO_ENABLE_CXX11 "-std=c++11" "CXX" "")

#============================================================
# Notes for bulding clang 3.5.0
#============================================================
# Attempt to install clang/llvm 3.5.0 on ccscs7
# - http://linuxdeveloper.blogspot.com/2012/12/building-llvm-32-from-source.html
# - In the source tree, create this symlink: llvm-3.5.0.src/tools/clang -> ../../cfe-3.5.0.src
# - Used modules for cmake, gcc/4.8.1.
# - ../llvm/configure \
#   --prefix=/ccs/codes/radtran/vendors/llvm-3.5.0 --enable-shared \
#   --with-gcc-toolchain=/ccs/codes/radtran/vendors/gcc-4.8.1 \
# - make -j 8 [install]

#------------------------------------------------------------------------------#
# End config/apple-clang.cmake
#------------------------------------------------------------------------------#
