#-----------------------------*-cmake-*----------------------------------------#
# file   config/unix-intel.cmake
# author Kelly Thompson
# date   2010 Nov 1
# brief  Establish flags for Linux64 - Intel C++
# note   Copyright (C) 2010-2013 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

# History
# ----------------------------------------
# 07/20/2011 - Use -O3 for Release builds but reduce -fp-model from
#              strict to precise to eliminate warning 1678.
# 11/18/2013 - For RELEASE builds, begin using -fp-model precise
#              -fp-speculation safe.  Jayenne sees about 10%
#              performance bump.

#
# Compiler Flags
#

if( NOT CXX_FLAGS_INITIALIZED )
   set( CXX_FLAGS_INITIALIZED "yes" CACHE INTERNAL "using draco settings." )

  set( CMAKE_C_FLAGS                "-w1 -vec-report0 -diag-disable remark -shared-intel -ftz" )
  set( CMAKE_C_FLAGS_DEBUG          "-g -O0 -inline-level=0 -ftrapuv -check=uninit -DDEBUG")
  set( CMAKE_C_FLAGS_RELEASE        "-O3 -ip -fp-speculation fast -fp-model fast -pthread -DNDEBUG" )
  set( CMAKE_C_FLAGS_MINSIZEREL     "${CMAKE_C_FLAGS_RELEASE}" )
  set( CMAKE_C_FLAGS_RELWITHDEBINFO "-g -debug inline-debug-info -O3 -ip -fp  -pthread -fp-model precise -fp-speculation safe" )

  set( CMAKE_CXX_FLAGS                "${CMAKE_C_FLAGS}" )
  set( CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_C_FLAGS_DEBUG} -early-template-check")
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

# Optional compiler flags
toggle_compiler_flag( HAVE_MIC                 "-mmic"        "C;CXX;EXE_LINKER" "")

# Options -mmic and -xHost are not compatible.  This check needs to
# appear after "-mmic" has been added to the compiler flags so that
# xhost will report that it is not available.
include(CheckCCompilerFlag)
check_c_compiler_flag(-xHost HAS_XHOST)

toggle_compiler_flag( HAS_XHOST                "-xHost"       "C;CXX"  "")
toggle_compiler_flag( USE_OPENMP               "-openmp"      "C;CXX;EXE_LINKER" "")
toggle_compiler_flag( DRACO_ENABLE_CXX11       "-std=c++0x"   "CXX" "")
toggle_compiler_flag( DRACO_ENABLE_C99         "-std=c99"     "C"   "")
toggle_compiler_flag( DRACO_ENABLE_STRICT_ANSI "-strict-ansi" "C;CXX" "")

#------------------------------------------------------------------------------#
# End config/unix-intel.cmake
#------------------------------------------------------------------------------#
