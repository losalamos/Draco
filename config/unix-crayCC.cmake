#-----------------------------*-cmake-*----------------------------------------#
# file   config/unix-crayCC.cmake
# author Kelly Thompson
# date   2010 Nov 1
# brief  Establish flags for Linux64 - Cray C/C++
# note   Copyright (C) 2016-2019 Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

# History
# ----------------------------------------
# 2015/09/28 - First cut.
# 6/13/2016  - IPO settings moved to compilerEnv.cmake
#              (CMAKE_INTERPROCEDURAL_OPTIMIZATION=ON).

#
# Compiler flag checks
#
if( CMAKE_CXX_COMPILER_VERSION LESS "8.4" )
  message( FATAL_ERROR "Cray C++ prior to 8.4 does not support C++11.
Try: module use --append ~mrberry/modulefiles
     module swap cce cce/8.4.0.233.
")
endif()

#
# Compiler Flags
#

# As of cmake-3.12.1, the cray flags for C++14 were not available in cmake.
if( NOT DEFINED CMAKE_CXX14_STANDARD_COMPILE_OPTION )
  set( CMAKE_CXX14_STANDARD_COMPILE_OPTION "-hstd=c++14" CACHE STRING "internal"
    FORCE)
endif()

if( NOT CXX_FLAGS_INITIALIZED )
   set( CXX_FLAGS_INITIALIZED "yes" CACHE INTERNAL "using draco settings." )

  set( CMAKE_C_FLAGS                "-DR123_USE_GNU_UINT128=0" )
  set( CMAKE_C_FLAGS_DEBUG          "-g -O0 -DDEBUG")
  #if( HAVE_MIC )
    # For floating point consistency with Xeon when using Intel 15.0.090 + Intel
    # MPI 5.0.2
    #set( CMAKE_C_FLAGS_DEBUG        "${CMAKE_C_FLAGS_DEBUG} -fp-model precise -fp-speculation safe" )
  #endif()
  set( CMAKE_C_FLAGS_RELEASE        "-O3 -DNDEBUG" )
  set( CMAKE_C_FLAGS_MINSIZEREL     "${CMAKE_C_FLAGS_RELEASE}" )
  set( CMAKE_C_FLAGS_RELWITHDEBINFO "-g -O3 -DNDEBUG" )

  set( CMAKE_CXX_FLAGS                "${CMAKE_C_FLAGS}")
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

# Optional compiler flags
# toggle_compiler_flag( HAVE_MIC                 "-mmic"        "C;CXX;EXE_LINKER" "")

# Options -mmic and -xHost are not compatible.  This check needs to
# appear after "-mmic" has been added to the compiler flags so that
# xhost will report that it is not available.
#include(CheckCCompilerFlag)
#check_c_compiler_flag(-xHost HAS_XHOST)
# If this is trinitite/trinity, do not use -xHost because front and back ends are different
# architectures. Instead use -xCORE-AVX2 (the default).
#if( NOT ${SITENAME} STREQUAL "Trinitite" )
#  toggle_compiler_flag( HAS_XHOST "-xHost" "C;CXX"  "")
#endif()
#toggle_compiler_flag( OPENMP_FOUND ${OpenMP_C_FLAGS} "C;CXX;EXE_LINKER" "" )
toggle_compiler_flag( OPENMP_FOUND ${OpenMP_C_FLAGS} "C;CXX" "" )

#------------------------------------------------------------------------------#
# End config/unix-crayCC.cmake
#------------------------------------------------------------------------------#
