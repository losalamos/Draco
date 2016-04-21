#-----------------------------*-cmake-*----------------------------------------#
# file   config/unix-intel.cmake
# author Kelly Thompson
# date   2010 Nov 1
# brief  Establish flags for Linux64 - Intel C++
# note   Copyright (C) 2016 Los Alamos National Security, LLC.
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
# Compiler flag checks
#
include(platform_checks)
query_openmp_availability()

#
# Compiler Flags
#

if( NOT CXX_FLAGS_INITIALIZED )
   set( CXX_FLAGS_INITIALIZED "yes" CACHE INTERNAL "using draco settings." )

  # [KT 2015-07-10] I would like to turn on -w2, but this generates
  #    many warnings from Trilinos headers that I can't suppress easily
  #    (warning 191: type qualifier is meaningless on cast type)
  # [KT 2015-07-10] -diag-disable 11060 -- disable warning that is
  #    issued when '-ip' is turned on and a library has no symbols (this
  #    occurs when capsaicin links some trilinos libraries.)
  set( CMAKE_C_FLAGS                "-w1 -vec-report0 -diag-disable remark -shared-intel -ftz -diag-disable 11060" )
  set( CMAKE_C_FLAGS_DEBUG          "-g -O0 -inline-level=0 -ftrapuv -check=uninit -DDEBUG")
  if( HAVE_MIC )
    # For floating point consistency with Xeon when using Intel 15.0.090 + Intel MPI 5.0.2
    set( CMAKE_C_FLAGS_DEBUG        "${CMAKE_C_FLAGS_DEBUG} -fp-model precise -fp-speculation safe" )
  endif()
  set( CMAKE_C_FLAGS_RELEASE        "-O3 -ip -fp-speculation fast -fp-model fast -pthread -DNDEBUG" )
  set( CMAKE_C_FLAGS_MINSIZEREL     "${CMAKE_C_FLAGS_RELEASE}" )
  set( CMAKE_C_FLAGS_RELWITHDEBINFO "-g -debug inline-debug-info -O3 -ip -fp  -pthread -fp-model precise -fp-speculation safe" )

  set( CMAKE_CXX_FLAGS                "${CMAKE_C_FLAGS} -std=c++11" )
  set( CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_C_FLAGS_DEBUG} -early-template-check")
  set( CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_C_FLAGS_RELEASE}")
  set( CMAKE_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_RELEASE}")
  set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}" )

   # Use C99 standard.
   set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")

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
# If this is trinitite/trinity, do not use -xHost because front and back ends are different
# architectures. Instead use -xCORE-AVX2 (the default).
if( NOT ${SITENAME} STREQUAL "Trinitite" )
  toggle_compiler_flag( HAS_XHOST "-xHost" "C;CXX"  "")
endif()
#toggle_compiler_flag( OPENMP_FOUND ${OpenMP_C_FLAGS} "C;CXX;EXE_LINKER" "" )
toggle_compiler_flag( OPENMP_FOUND ${OpenMP_C_FLAGS} "C;CXX" "" )

#
# Sanity checks
#

# On Moonlight, Intel-16 requires MKL-11.3
if( CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 15.1 AND DEFINED ENV{MKLROOT} )
  if( NOT $ENV{MKLROOT} MATCHES "2016" )
    message( FATAL_ERROR "Intel-16 requires MKL-11.3+.")
  endif()
endif()

#------------------------------------------------------------------------------#
# End config/unix-intel.cmake
#------------------------------------------------------------------------------#
