#-----------------------------*-cmake-*----------------------------------------#
# file   device/test/CPU_Device.cmake
# brief  Instructions for building device/test Makefile.
# note   Copyright (C) 2016-2019 Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_sources
  gpu_hello_rt_api.cu
  gpu_device_info.cu
  gpu_dual_call_test.cu
)

set(cuda_headers
  basic_kernels.hh
  Dual_Call.hh
)

set(cuda_sources
  basic_kernels.cu
  Dual_Call.cu
)

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

# Stuff cuda code into a test library.
add_component_library(
   TARGET       Lib_device_test
   TARGET_DEPS  Lib_dsxx Lib_device
   LIBRARY_NAME device_test
   LIBRARY_TYPE STATIC
   SOURCES      "${cuda_sources}"
   HEADERS      "${cuda_headers}" )
set_property(TARGET Lib_device_test PROPERTY CUDA_SEPARABLE_COMPILATION ON)

# We need to explicitly state that we need all CUDA files in the particle
# library to be built with -dc as the member functions could be called by
# other libraries and executables
set_target_properties( Lib_device_test PROPERTIES
 CUDA_SEPARABLE_COMPILATION ON )

#target_include_directories( Ut_gpu_hello_rt_api_exe
# PRIVATE $<BUILD_INTERFACE:${draco_src_dir_SOURCE_DIR}> )

#          $<BUILD_INTERFACE:${dsxx_BINARY_DIR}> )
#  PUBLIC $<BUILD_INTERFACE:${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}> )


#   ${CUDA_CUDA_LIBRARY} )
#   ${CUDA_TOOLKIT_INCLUDE}

# ---------------------------------------------------------------------------- #
# Register Unit tests
# ---------------------------------------------------------------------------- #

 set( test_deps
    Lib_device
    Lib_device_test)

 add_scalar_tests(
    SOURCES  "${test_sources}"
    DEPS     "${test_deps}"
    )

#add_test( NAME device_gpu_hello_rt_api
#  COMMAND $<TARGET_FILE:Ut_gpu_hello_rt_api_exe> )
set_tests_properties( device_gpu_hello_rt_api PROPERTIES
   PASS_REGULAR_EXPRESSION ".*[Tt]est: PASSED"
   FAIL_REGULAR_EXPRESSION ".*[Tt]est: FAILED" )

# ---------------------------------------------------------------------------- #
# End device/test/GPU_Device.cmake
# ---------------------------------------------------------------------------- #
