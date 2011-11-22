#-----------------------------*-cmake-*----------------------------------------#
# file   config/unix-pgi.cmake
# author Kelly Thompson 
# date   2010 Nov 1
# brief  Establish flags for Linux64 - Intel C++
# note   Copyright © 2010 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

# NOTE: You may need to set TMPDIR to a location that the compiler can
#       write temporary files to.  Sometimes the default space on HPC
#       is not sufficient.

#
# Sanity Checks
# 

if( BUILD_SHARED_LIBS )
  message( FATAL_ERROR "Feature not available - yell at KT." )
endif( BUILD_SHARED_LIBS )

# Cannot use strict ansi flags on RedStorm
option( ENABLE_STRICT_ANSI "Turn on strict ANSI compliance?" ON )
if( ${ENABLE_STRICT_ANSI} )
   set( STRICT_ANSI_FLAGS "-Xa -A --no_using_std" )
endif()


#
# C++ libraries required by Fortran linker
# 

#
# config.h settings
#

execute_process(
  COMMAND ${CMAKE_C_COMPILER} -V
  OUTPUT_VARIABLE ABS_C_COMPILER_VER
  )
string( REGEX REPLACE "Copyright.*" " " 
  ABS_C_COMPILER_VER ${ABS_C_COMPILER_VER} )
string( STRIP ${ABS_C_COMPILER_VER} ABS_C_COMPILER_VER )

execute_process(
  COMMAND ${CMAKE_CXX_COMPILER} -V
  OUTPUT_VARIABLE ABS_CXX_COMPILER_VER
  )
string( REGEX REPLACE "Copyright.*" " " 
  ABS_CXX_COMPILER_VER ${ABS_CXX_COMPILER_VER} )
string( STRIP ${ABS_CXX_COMPILER_VER} ABS_CXX_COMPILER_VER )

#
# Compiler Flags
# 

# http://www.pgroup.com/support/compile.htm

# Flags from Draco autoconf build system:
# -Xa
# -A                       ansi
# --no_using_std           Enable (disable) implicit use of the std
#                          namespace when standard header files are
#                          included. 
# --diag_suppress 940      Suppress warning #940
# --diag_suppress 11           "      "     # 11
# -DNO_PGI_OFFSET
# -Kieee                   Perform floating-point operations in strict
#                          conformance with the IEEE 754 standard. 
# --no_implicit_include    Disable implicit inclusion of source files
#                          as a method of finding definitions of
#                          template entities to be instantiated. 
# -Mdaz                    Enable (disable) mode to treat denormalized
#                          floating point numbers as zero.  -Mdaz is
#                          default for -tp p7-64 targets; -Mnodaz is
#                          default otherwise. 
# -pgf90libs               Link-time option to add the pgf90 runtime
#                          libraries, allowing mixed-language programming. 
# -Mipa                    Enable and specify options for
#                          InterProcedural Analysis (IPA). 
# -Mnoframe                Don't setup a true stack frame pointer for
#                          functions. allows slightly more efficient
#                          operation when a stack frame is not needed.
# -Mlre                    Enable loop-carried redundancy elimination.
# -Mautoinline=levels:n    Enable inlining of functions with the
#                          inline attribute up to n levels deep.  The
#                          default is to inline up to 5 levels. 
# -Mvect=sse               Use SSE, SSE2, 3Dnow, and prefetch
#                          instructions in loops where possible. 
# -Mcache_align            Align unconstrained data objects of size
#                          greater than or equal to 16 bytes on
#                          cache-line boundaries.  An unconstrained
#                          object is a variable or array that is not a
#                          member of an aggregate structure or common
#                          block, is not allocatable, and is not an
#                          automatic array.
# -Mflushz                 Set SSE to flush-to-zero mode.


if( NOT CXX_FLAGS_INITIALIZED )
   set( CXX_FLAGS_INITIALIZED "yes" CACHE INTERNAL "using draco settings." )

   set( CMAKE_C_FLAGS                "-Kieee -Mdaz -pgf90libs" )
   set( CMAKE_C_FLAGS_DEBUG          "-g -O0") # -DDEBUG") 
   set( CMAKE_C_FLAGS_RELEASE        "-O3 -DNDEBUG" ) # -O4
   set( CMAKE_C_FLAGS_MINSIZEREL     "${CMAKE_C_FLAGS_RELEASE}" )
   set( CMAKE_C_FLAGS_RELWITHDEBINFO "-O3 -DNDEBUG -gopt" )

   set( CMAKE_CXX_FLAGS                "${CMAKE_C_FLAGS} ${STRICT_ANSI_FLAGS} --no_implicit_include --diag_suppress 940 --diag_suppress 11 --diag_suppress 450 -DNO_PGI_OFFSET" )
   set( CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_C_FLAGS_DEBUG}")
   set( CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_C_FLAGS_RELEASE} -Munroll=c:10 -Mautoinline=levels:10 -Mvect=sse -Mflushz -Mlre")

# -Mipa=fast,inline
# -zc_eh 
# -Msmartalloc
# -tp x64      Create a PGI Unified Binary which functions correctly
#              on and is optimized for both Intel and AMD processors. 
# -Mprefetch   Control generation of prefetch instructions to improve
#              memory performance in compute-intensive loops. 

# -Mnoframe (we use this to debug crashed programs).
# -Mcache_align (breaks some tests in wedgehog)
# -Msafeptr (breaks some array operations in MatRA).
   set( CMAKE_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_RELEASE}")
   set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}" )
ENDIF()

string( TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UPPER )

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
# End config/unix-pgi.cmake
#------------------------------------------------------------------------------#
