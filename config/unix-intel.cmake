#-----------------------------*-cmake-*----------------------------------------#
# file   config/unix-intel.cmake
# author Kelly Thompson
# date   2010 Nov 1
# brief  Establish flags for Linux64 - Intel C++
# note   Copyright (C) 2016-2019 Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

# History
# ----------------------------------------
# 07/20/2011 - Use -O3 for Release builds but reduce -fp-model from strict to
#              precise to eliminate warning 1678.
# 11/18/2013 - For RELEASE builds, begin using -fp-model precise -fp-speculation
#              safe.  Jayenne sees about 10% performance bump.
# 06/13/2016 - IPO settings moved to compilerEnv.cmake
#              (CMAKE_INTERPROCEDURAL_OPTIMIZATION=ON).

#
# Compiler Flags
#

if( NOT CXX_FLAGS_INITIALIZED )
   set( CXX_FLAGS_INITIALIZED "yes" CACHE INTERNAL "using draco settings." )

  # [KT 2015-07-10] I would like to turn on -w2, but this generates many
  #    warnings from Trilinos headers that I can't suppress easily (warning 191:
  #    type qualifier is meaningless on cast type)
  # [KT 2015-07-10] -diag-disable 11060 -- disable warning that is issued when
  #    '-ip' is turned on and a library has no symbols (this occurs when
  #    capsaicin links some trilinos libraries.)
  set( CMAKE_C_FLAGS
    "-w1 -vec-report0 -diag-disable remark -shared-intel -no-ftz -fma -diag-disable 11060" )
  set( CMAKE_C_FLAGS_DEBUG
    "-g -O0 -inline-level=0 -ftrapuv -check=uninit -fp-model precise -fp-speculation safe -debug inline-debug-info -fno-omit-frame-pointer -DDEBUG")
  set( CMAKE_C_FLAGS_RELEASE
    "-O3 -fp-speculation fast -fp-model precise -pthread -DNDEBUG" )
  # [KT 2017-01-19] On KNL, -fp-model fast changes behavior significantly for
  # IMC. Revert to -fp-model precise.
  if( "$ENV{CRAY_CPU_TARGET}" STREQUAL "mic-knl" )
    string( REPLACE "-fp-model fast" "-fp-model precise" CMAKE_C_FLAGS_RELEASE
      ${CMAKE_C_FLAGS_RELEASE} )
  endif()

  set( CMAKE_C_FLAGS_MINSIZEREL
    "${CMAKE_C_FLAGS_RELEASE}" )
  set( CMAKE_C_FLAGS_RELWITHDEBINFO
    "-g -debug inline-debug-info -O3 -pthread -fp-model precise -fp-speculation safe -fno-omit-frame-pointer" )

  set( CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS}")
  set( CMAKE_CXX_FLAGS_DEBUG
    "${CMAKE_C_FLAGS_DEBUG} -early-template-check")
  set( CMAKE_CXX_FLAGS_RELEASE
    "${CMAKE_C_FLAGS_RELEASE}")
  set( CMAKE_CXX_FLAGS_MINSIZEREL
    "${CMAKE_CXX_FLAGS_RELEASE}")
  set( CMAKE_CXX_FLAGS_RELWITHDEBINFO
    "${CMAKE_C_FLAGS_RELWITHDEBINFO}" )

   # Use C99 standard.
   set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")

endif()

find_library( INTEL_LIBM m )
mark_as_advanced( INTEL_LIBM )

##---------------------------------------------------------------------------##
# Ensure cache values always match current selection
##---------------------------------------------------------------------------##
set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS}" CACHE STRING "compiler flags" FORCE )
set( CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}" CACHE STRING "compiler flags"
  FORCE )
set( CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}" CACHE STRING
  "compiler flags" FORCE )
set( CMAKE_C_FLAGS_MINSIZEREL "${CMAKE_C_FLAGS_MINSIZEREL}" CACHE STRING
  "compiler flags" FORCE )
set( CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}" CACHE STRING
  "compiler flags" FORCE )

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" CACHE STRING "compiler flags" FORCE )
set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}" CACHE STRING
  "compiler flags" FORCE )
set( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING
  "compiler flags" FORCE )
set( CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL}" CACHE STRING
  "compiler flags" FORCE )
set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}" CACHE
  STRING "compiler flags" FORCE )

# If this is a Cray, the compile wrappers take care of any xHost flags that are
# needed.
if( NOT CRAY_PE )
# include(CheckCCompilerFlag)
# check_c_compiler_flag(-xHost HAS_XHOST)
 set( HAS_XHOST TRUE )
 toggle_compiler_flag( HAS_XHOST "-xHost" "C;CXX" "")
#else()
 # -craype-verbose
endif()
toggle_compiler_flag( OPENMP_FOUND ${OpenMP_C_FLAGS} "C;CXX" "" )

#------------------------------------------------------------------------------#
# End config/unix-intel.cmake
#------------------------------------------------------------------------------#
