#-----------------------------*-cmake-*----------------------------------------#
# file   config/unix-intel.cmake
# author Kelly Thompson 
# date   2010 Nov 1
# brief  Establish flags for Linux64 - Intel C++
# note   Copyright © 2010 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

# History
# ----------------------------------------
# 7/20/11 - Use -O3 for Release builds but reduce -fp-model from
#           strict to precise to eliminate warning 1678.



#
# Sanity Checks
# 

if( NOT __LINUX_COMPILER_INTEL )
  message( FATAL_ERROR "If CXX is not Intel C++, then we shouldn't have ended up here.  Something is really wrong with the build system. " )
endif()

#
# Compiler Flags
# 

# Try 'icpc -help':
# -inline-level=<n> control inline expansion (same as -Ob<n>)
#    n=0  disables inlining
#    n=1  inline functions declared with __inline, and perform C++ inlining
#    n=2  inline any function, at the compiler's discretion (same as -ip)
# -O3    enable -O2 plus more aggressive optimizations that may not improve
#        performance for all programs
# -O0    disable optimizations
# -g     Include debug information
# -ip    enable single-file IP optimizations (within files)
# -ipo   enable multi-file IP optimizations (within files)
# -ansi  equivalent to GNU -ansi
# -fp-model <name>    enable <name> floating point model variation
#            [no-]except - enable/disable floating point semantics
#            double      - rounds intermediates in 53-bit (double) precision
#            extended    - rounds intermediates in 64-bit (extended) precision
#            fast[=1|2]  - enables more aggressive floating point optimizations
#            precise     - allows value-safe optimizations
#            source      - enables intermediates in source precision
#            strict      - enables -fp-model precise -fp-model except, disables
#                          contractions and enables pragma stdc fenv_access
# -w<n>      control diagnostics:
#            n=0 display errors (same as -w)
#            n=1 display warnings and errors (DEFAULT)
#            n=2 display remarks, warnings, and errors
# -shared-intel Causes Intel-provided libraries to be linked in
#            dynamically.  The should eliminate the need to link
#            against libm for every library.

# Suppressions 
# -wd<L1>[,<L2>,...] disable diagnostics L1 through LN
# Warning #1678: cannot enable speculation unless fenv_access and 
#                exception_semantics are disabled

if( NOT CXX_FLAGS_INITIALIZED )
   set( CXX_FLAGS_INITIALIZED "yes" CACHE INTERNAL "using draco settings." )

  set( CMAKE_C_FLAGS                "-fpic -w1 -vec-report0 -diag-disable remark -shared-intel" )
  set( CMAKE_C_FLAGS_DEBUG          "-g -O0 -inline-level=0 -ftrapuv -DDEBUG") 
  set( CMAKE_C_FLAGS_RELEASE        "-O3 -inline-level=1 -ip -fpe0 -fp-model precise -fp-speculation strict -ftz -pthread -DNDEBUG" )
  set( CMAKE_C_FLAGS_MINSIZEREL     "${CMAKE_C_FLAGS_RELEASE}" )
  set( CMAKE_C_FLAGS_RELWITHDEBINFO "-O3 -inline-level=0 -ip -DNDEBUG" )

  set( CMAKE_CXX_FLAGS                "${CMAKE_C_FLAGS}" )
  set( CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_C_FLAGS_DEBUG} -strict-ansi")
  set( CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_C_FLAGS_RELEASE}")
  set( CMAKE_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_RELEASE}")
  set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}" )
endif()

find_library( INTEL_LIBM m )
mark_as_advanced( INTEL_LIBM )


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
