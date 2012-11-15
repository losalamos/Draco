#-----------------------------*-cmake-*----------------------------------------#
# file   config/unix-xl.cmake
# author Gabriel Rockefeller 
# date   2012 Nov 1
# brief  Establish flags for Linux64 - IBM XL C++
# note   Copyright (C) 2012 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

# Notable xlc/xlc++ options:

# -qinfo=all      Enable all warnings (and informational messages).
# -qflag=i:w      Send informational-level messages to a listing file,
#                 if one is requested, but only send warning-level or
#                 more severe messages (i.e., errors) to the terminal.
# -qcheck         Enable runtime bounds, null-pointer, and
#                 divide-by-zero checks.
# -O3             -O2 plus memory- and compile-time-intensive
#                 operations that can alter the semantics of programs.
# -qstrict=nans:operationprecision
#                 Disable optimizations at -O3 and above that may
#                 produce incorrect results in the presence of NaNs,
#                 or that produce approximate results for individual
#                 floating-point operations.
# -qhot=novector  Enable high-order transformations during
#                 optimization.  LLNL recommends novector to disable
#                 gathering math intrinsics into separate vector math
#                 library calls (because it's typically better to let
#                 those instructions intermix with other
#                 floating-point operations, when using SIMD
#                 instructions).
# -qsimd=auto     Enable automatic generation of SIMD instructions, to
#                 take advantage of BG/Q-specific Quad Processing
#                 eXtension (QPX) units.
# -qnostaticlink  Allow linking with shared libraries.
# -qsmp=omp       Enable parallelization using OpenMP pragmas.

# Suppressions:
# -qsuppress  Suppress specific informational or warning messages.
#  1540-0072: The attribute [...] is not supported on the target
#             platform. The attribute is ignored.

if( NOT CXX_FLAGS_INITIALIZED )
   set( CXX_FLAGS_INITIALIZED "yes" CACHE INTERNAL "using draco settings." )

   set( CMAKE_C_FLAGS                "-qinfo=all -qflag=i:w -qsuppress=1540-0072" )
   set( CMAKE_C_FLAGS_DEBUG          "-g -O0 -qcheck -DDEBUG")
   set( CMAKE_C_FLAGS_RELEASE        "-O3 -qhot=novector -qsimd=auto -qstrict=nans:operationprecision -DNDEBUG" )
   set( CMAKE_C_FLAGS_MINSIZEREL     "${CMAKE_C_FLAGS_RELEASE}" )
   set( CMAKE_C_FLAGS_RELWITHDEBINFO "-g -O3 -qhot=novector -qsimd=auto -qstrict=nans:operationprecision" )

   set( CMAKE_CXX_FLAGS                "${CMAKE_C_FLAGS}" )
   set( CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_C_FLAGS_DEBUG}")
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

toggle_compiler_flag( DRACO_SHARED_LIBS "-qnostaticlink" "EXE_LINKER" "")
toggle_compiler_flag( USE_OPENMP       "-qsmp=omp" "C;CXX;EXE_LINKER" "")
toggle_compiler_flag( DRACO_ENABLE_CXX11 "-qlanglvl=extended0x" "CXX" "")
toggle_compiler_flag( DRACO_ENABLE_STRICT_ANSI "-qlanglvl=stdc89" "C" "")

#------------------------------------------------------------------------------#
# End config/unix-xl.cmake
#------------------------------------------------------------------------------#
