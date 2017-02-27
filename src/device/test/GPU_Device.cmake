#-----------------------------*-cmake-*----------------------------------------#
# file   device/test/CMakeLists.txt
# brief  Instructions for building device/test Makefile.
# note   Copyright (C) 2016-2017 Los Alamos National Security
#        All rights reserved.
# version $Id$
#------------------------------------------------------------------------------#

if( NOT EXISTS ${CUDA_NVCC_EXECUTABLE} )
   find_package(CUDA)
endif()

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( test_sources
   # gpu_hello_rt_api.cc
   gpu_hello_driver_api.cc
   )
set( cuda_sources
   gpu_kernel.cu
   vector_add.cu
   )

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories(
   ${PROJECT_SOURCE_DIR}        # headers for tests
   ${PROJECT_SOURCE_DIR}/..     # headers for package
   ${PROJECT_BINARY_DIR}/..     # config.h
   ${CUDA_TOOLKIT_INCLUDE}
   )

# ---------------------------------------------------------------------------- #
# Compile .cu files into .ptx files.
# ---------------------------------------------------------------------------- #

set( CUDA_BUILD_CUBIN ON )
cuda_compile( generated_files ${cuda_sources} )

# move and rename the cubin.txt output file from the CMakeLists
# directory to the CWD.
# CMakeFiles/cuda_compile.dir/cuda_compile_generated_gpu_kernel.cu.o.cubin.txt
#    -> gpu_kernel.cubin
unset( cubinfiles_list )
foreach( cubinfile ${generated_files} )
   get_filename_component( sh_cubinfile ${cubinfile} NAME_WE )
   string( REPLACE "cuda_compile_generated_" "" sh_cubinfile
      ${sh_cubinfile} )
   set( short_cubinfile "${sh_cubinfile}.cubin" )
   add_custom_command( OUTPUT ${short_cubinfile}
      COMMAND ${CMAKE_COMMAND} -E rename ${cubinfile}.cubin.txt ${short_cubinfile}
      DEPENDS ${cubinfile} )
   list(  APPEND cubinfiles_list ${short_cubinfile} )
endforeach()
add_custom_target( device_test_build_cubin_files ALL
   DEPENDS ${cubinfiles_list} )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

# Binary using Runtime API (no link to libraries other than libcuda)
cuda_add_executable( gpu_hello_rt_api_exe
   ${PROJECT_SOURCE_DIR}/gpu_hello_rt_api.cc )
target_link_libraries( gpu_hello_rt_api_exe Lib_dsxx )

# Binary using Device API and links to libcudahello.
#   add_executable( gpu_hello_driver_api_exe
#      ${PROJECT_SOURCE_DIR}/gpu_hello_driver_api.cc
#      )
#   target_link_libraries( gpu_hello_driver_api_exe
#      Lib_device
#      Lib_dsxx
#      ${CUDA_CUDA_LIBRARY}
#      )

# ---------------------------------------------------------------------------- #
# Register Unit tests
# ---------------------------------------------------------------------------- #

set( test_deps
   Lib_device
   ${CUDA_CUDA_LIBRARY} )
add_scalar_tests(
   SOURCES  "${test_sources}"
   DEPS     "${test_deps}"
   )

add_test(
   NAME    device_gpu_hello_rt_api_exe
   COMMAND $<TARGET_FILE:gpu_hello_rt_api_exe>
   )
set_tests_properties( device_gpu_hello_rt_api_exe
   PROPERTIES
   PASS_REGULAR_EXPRESSION ".*[Tt]est: PASSED"
   FAIL_REGULAR_EXPRESSION ".*[Tt]est: FAILED"
   )

set( extra_clean_files
   cuda_compile_generated_gpu_kernel.ptx
   cuda_compile_generated_vector_add.ptx
   gpu_kernel.cubin
   gpu_kernel.ptx
   vector_add.cubin
   vector_add.ptx
   )
set_directory_properties(
   PROPERTIES
   ADDITIONAL_MAKE_CLEAN_FILES "${extra_clean_files}" )
