#-----------------------------*-cmake-*----------------------------------------#
# file   config/unix-cuda.cmake
# brief  Establish flags for Unix/Linux - Cuda
# note   Copyright (C) 2020 Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

include_guard(GLOBAL)

# Note: In config/component_macros.cmake, the build system sets flags for
# 1) the language standard (C++14, C99, etc)
# 2) interprocedural optimization.

# Notes:
# ----------------------------------------
# Useful options that could be added to aid debugging
# - https://devblogs.nvidia.com/tag/cuda/
# - https://devblogs.nvidia.com/building-cuda-applications-cmake/

# Discover CUDA compute capabilities.
if( NOT DEFINED Draco_CUDA_ARCH )
  if( NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/config )
    file(MAKE_DIRECTORY  ${CMAKE_CURRENT_BINARY_DIR}/config )
  endif()
  set(OUTPUTFILE ${CMAKE_CURRENT_BINARY_DIR}/config/cuda_script)
  set(CUDAFILE ${CMAKE_CURRENT_SOURCE_DIR}/config/query_gpu.cu)
  execute_process( COMMAND ${CMAKE_CUDA_COMPILER}
    -lcuda ${CUDAFILE} -o ${OUTPUTFILE} )
  execute_process( COMMAND ${OUTPUTFILE}
    RESULT_VARIABLE CUDA_RETURN_CODE
    OUTPUT_VARIABLE Draco_CUDA_ARCH )
  if( NOT ${CUDA_RETURN_CODE} EQUAL 0 )
    message( FATAL_ERROR "Unable to determine target Cuda arch." )
  endif()
  unset(OUTPUTFILE)
  unset(CUDAFILE)
  unset(CUDA_RETURN_CODE)
endif()

#
# Compiler Flags
#
# -Xcudafe --diag_suppress=1427: suppress "warning: offsetof applied to non-POD"

if( NOT CUDA_FLAGS_INITIALIZED )
  set( CUDA_FLAGS_INITIALIZED "yes" CACHE INTERNAL "using draco settings." )

  set( CMAKE_CUDA_FLAGS                "${Draco_CUDA_ARCH} -g -G --expt-relaxed-constexpr")
  set( CMAKE_CUDA_FLAGS_DEBUG          "-O0 -Xcudafe --display_error_number -Xcudafe --diag_suppress=1427")
  set( CMAKE_CUDA_FLAGS_RELEASE        "-O2")
  set( CMAKE_CUDA_FLAGS_MINSIZEREL     "-O2")
  set( CMAKE_CUDA_FLAGS_RELWITHDEBINFO "-O2 --generate-line-info")

endif()

##---------------------------------------------------------------------------##
# Ensure cache values always match current selection
##---------------------------------------------------------------------------##
set( CMAKE_CUDA_FLAGS                "${CMAKE_CUDA_FLAGS}"                CACHE
     STRING "compiler flags" FORCE )
set( CMAKE_CUDA_FLAGS_DEBUG          "${CMAKE_CUDA_FLAGS_DEBUG}"          CACHE
     STRING "compiler flags" FORCE )
set( CMAKE_CUDA_FLAGS_RELEASE        "${CMAKE_CUDA_FLAGS_RELEASE}"        CACHE
     STRING "compiler flags" FORCE )
set( CMAKE_CUDA_FLAGS_MINSIZEREL     "${CMAKE_CUDA_FLAGS_MINSIZEREL}"     CACHE
     STRING "compiler flags" FORCE )
set( CMAKE_CUDA_FLAGS_RELWITHDEBINFO "${CMAKE_CUDA_FLAGS_RELWITHDEBINFO}" CACHE
     STRING "compiler flags" FORCE )

#------------------------------------------------------------------------------#
# End config/unix-cuda.cmake
#------------------------------------------------------------------------------#
