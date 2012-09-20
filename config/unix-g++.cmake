#-----------------------------*-cmake-*----------------------------------------#
# file   config/unix-g++.cmake
# brief  Establish flags for Windows - MSVC
# note   Copyright (C) 2010-2012 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

# Useful options that could be added to aid debugging
# http://gcc.gnu.org/wiki/Atomic/GCCMM
# http://gcc.gnu.org/wiki/TransactionalMemory
# http://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html


#
# Sanity Checks
# 

if( NOT CMAKE_COMPILER_IS_GNUCC )
  message( FATAL_ERROR "If CC is not GNUCC, then we shouldn't have ended up here.  Something is really wrong with the build system. " )
endif( NOT CMAKE_COMPILER_IS_GNUCC )

#
# Compiler flag check
#
include(CheckCCompilerFlag)
check_c_compiler_flag(-march=native HAS_MARCH_NATIVE)

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
set( DBS_C_COMPILER_VER "gcc ${CMAKE_CXX_COMPILER_VERSION}" )
set( DBS_CXX_COMPILER_VER "g++ ${CMAKE_CXX_COMPILER_VERSION}" )

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

set( DRACO_ENABLE_STRICT_ANSI ON CACHE INTERNAL 
   "use strict ANSI flags, C98" FORCE )

#
# Compiler Flags
# 

# cmake-2.8.9+ -fPIC can be established by using
# set(CMAKE_POSITION_INDEPENDENT_CODE ON) 

# Flags from Draco autoconf build system:

if( NOT CXX_FLAGS_INITIALIZED )
   set( CXX_FLAGS_INITIALIZED "yes" CACHE INTERNAL "using draco settings." )

   set( CMAKE_C_FLAGS                "-Wcast-align -Wpointer-arith -Wall" )
   if (NOT ${CMAKE_GENERATOR} MATCHES Xcode AND HAS_MARCH_NATIVE)
      set( CMAKE_C_FLAGS                "${CMAKE_C_FLAGS} -march=native" )
   endif()
   set( CMAKE_C_FLAGS_DEBUG          "-g -fno-inline -fno-eliminate-unused-debug-types -O0 -Wextra -DDEBUG")
   set( CMAKE_C_FLAGS_RELEASE        "-O3 -funroll-loops -DNDEBUG" )
   set( CMAKE_C_FLAGS_MINSIZEREL     "${CMAKE_C_FLAGS_RELEASE}" )
   set( CMAKE_C_FLAGS_RELWITHDEBINFO "-O3 -g -fno-eliminate-unused-debug-types -Wextra -funroll-loops" )

   set( CMAKE_CXX_FLAGS                "${CMAKE_C_FLAGS}" ) 
   set( CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_C_FLAGS_DEBUG} -Woverloaded-virtual -Wno-long-long")
   set( CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_C_FLAGS_RELEASE}")
   set( CMAKE_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_RELEASE}")
   set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}" )

endif()

# Extra flags for gcc-4.6.2+
# -Wsuggest-attribute=[const|pure|noreturn]
if( "${DBS_CXX_COMPILER_VER_MAJOR}" GREATER 3  AND NOT  ${CMAKE_GENERATOR} MATCHES Xcode ) # 4
   if( "${DBS_CXX_COMPILER_VER_MINOR}" GREATER 5 ) # 4.6
      # include(CheckCXXCompilerFlag)
      # CHECK_CXX_COMPILER_FLAG( "-Wnoexcept, HAS_WNOEXCEPT)
      if( NOT "${CMAKE_CXX_FLAGS_DEBUG}" MATCHES "-Wsuggest-attribute=const" )
         set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wsuggest-attribute=const" )
      endif()
      if( NOT "${CMAKE_CXX_FLAGS_DEBUG}" MATCHES "-Wnoexcept" )
         set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wnoexcept" )
      endif()
   endif()
   if( "${DBS_CXX_COMPILER_VER_MINOR}" GREATER 6 ) # 4.7 
      # http://gcc.gnu.org/gcc-4.7/changes.html
      if( NOT "${CMAKE_CXX_FLAGS_DEBUG}" MATCHES "-Wunused-local-typedefs" )
         set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wunused-local-typedefs" )
      endif()
      # if( NOT "${CMAKE_CXX_FLAGS_DEBUG}" MATCHES "-Wzero-as-null-pointer-constant" )
      #    set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wzero-as-null-pointer-constant" )
      # endif()
   endif()
endif()

##---------------------------------------------------------------------------##

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

# Toggle for OpenMP
toggle_compiler_flag( USE_OPENMP         "-fopenmp"   "C;CXX;EXE_LINKER" )
# Toggle for C++11 support
# can use -std=c++11 with version 4.7+
if( CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.7 OR CMAKE_CXX_COMPILER VERSION_EQUAL 4.7)
   toggle_compiler_flag( DRACO_ENABLE_CXX11 "-std=c++11" "CXX") 
   set( DRACO_ENABLE_STRICT_ANSI OFF CACHE INTERNAL 
      "disable strict ANSI" FORCE)
elseif(  CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.3 OR CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL 4.3)
   toggle_compiler_flag( DRACO_ENABLE_CXX11 "-std=c++0x" "CXX") 
   set( DRACO_ENABLE_STRICT_ANSI OFF CACHE INTERNAL 
      "disable strict ANSI" FORCE)
else()
   if( DRACO_ENABLE_CXX11 )
      message(FATAL_ERROR "
C++11 requested via DRACO_ENABLE_CXX11=${DRACO_ENABLE_CXX11}.  
Therefore a gcc compiler with a version higher than 4.3 is needed.
Found gcc version ${GCC_VERSION}") 
   endif()
endif()

# Do we add '-ansi -pedantic'?
toggle_compiler_flag( DRACO_ENABLE_STRICT_ANSI "-ansi -pedantic" "CXX" )

# Notes for building gcc-4.7.1
# ../gcc-4.7.1/configure \
# --prefix=/ccs/codes/radtran/vendors/gcc-4.7.1 \
# --enable-language=c++,fortran,lto \
# --with-gmp-lib=/ccs/codes/radtran/vendors/Linux64/gmp-4.3.2/lib \
# --with-gmp-include=/ccs/codes/radtran/vendors/Linux64/gmp-4.3.2/include \
# --with-mpfr-lib=/ccs/codes/radtran/vendors/Linux64/mpfr-3.0.0/lib \
# --with-mpfr-include=/ccs/codes/radtran/vendors/Linux64/mpfr-3.0.0/include \
# --with-mpc-lib=/ccs/codes/radtran/vendors/Linux64/mpc-0.8.2/lib \
# --with-mpc-include=/ccs/codes/radtran/vendors/Linux64/mpc-0.8.2/include 
# --disable-multilib

# gmp-5.0.5
# ./configure --prefix=/usr/projects/draco/vendors/gmp-5.0.5
# make; make check; make install

# mpfr-3.1.1
# ./configure --prefix=/usr/projects/draco/vendors/mpfr-3.1.1 --with-gmp=/usr/projects/draco/vendors/gmp-5.0.5
# make; make check; make install

# mpc-0.8.2
# ./configure --prefix=/usr/projects/draco/vendors/mpc-0.8.2 --with-gmp=/usr/projects/draco/vendors/gmp-5.0.5 --with-mpfr=/usr/projects/draco/vendors/mpfr-3.1.1
# make; make check; make install

# Note: On some redhat systems you may need to hide (move) all shared
# libraries for gmp, mpfr and mpc before configuring/making gcc.

#------------------------------------------------------------------------------#
# End config/unix-g++.cmake
#------------------------------------------------------------------------------#
