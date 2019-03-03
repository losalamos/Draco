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
  # The low level Driver-API is not really supported by NVidia anymore.
  # gpu_hello_driver_api.cc
  )

# set( cuda_sources
#   gpu_kernel.cu
#   vector_add.cu
#   )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

# Stuff cuda code into a test library.
# add_library( Lib_device_test SHARED ${cuda_sources} )

#add_executable( Ut_gpu_hello_driver_api_exe gpu_hello_driver_api.cc )
#target_link_libraries( Ut_gpu_hello_driver_api_exe Lib_device )

add_executable( Ut_gpu_hello_rt_api_exe gpu_hello_rt_api.cu )
target_link_libraries( Ut_gpu_hello_rt_api_exe Lib_dsxx )
# target_link_libraries( Ut_gpu_hello_rt_api_exe Lib_device_test Lib_dsxx )
#target_link_libraries( Ut_gpu_hello_rt_api_exe Lib_device )

# We need to explicitly state that we need all CUDA files in the particle
# library to be built with -dc as the member functions could be called by
# other libraries and executables
set_target_properties( Ut_gpu_hello_rt_api_exe PROPERTIES
  CUDA_SEPARABLE_COMPILATION ON )
target_include_directories( Ut_gpu_hello_rt_api_exe
  PRIVATE $<BUILD_INTERFACE:${draco_src_dir_SOURCE_DIR}> )

#          $<BUILD_INTERFACE:${dsxx_BINARY_DIR}> )
#  PUBLIC $<BUILD_INTERFACE:${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}> )


#   ${CUDA_CUDA_LIBRARY} )
#   ${CUDA_TOOLKIT_INCLUDE}

# ---------------------------------------------------------------------------- #
# Register Unit tests
# ---------------------------------------------------------------------------- #

# set( test_deps
#    Lib_device
#    ${CUDA_CUDA_LIBRARY} )
# add_scalar_tests(
#    SOURCES  "${test_sources}"
#    DEPS     "${test_deps}"
#    )

add_test( NAME device_gpu_hello_rt_api
  COMMAND $<TARGET_FILE:Ut_gpu_hello_rt_api_exe> )
set_tests_properties( device_gpu_hello_rt_api PROPERTIES
   PASS_REGULAR_EXPRESSION ".*[Tt]est: PASSED"
   FAIL_REGULAR_EXPRESSION ".*[Tt]est: FAILED" )

# ---------------------------------------------------------------------------- #
# End device/test/GPU_Device.cmake
# ---------------------------------------------------------------------------- #
