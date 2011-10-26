#-----------------------------*-cmake-*----------------------------------------#
# file   device/test/CMakeLists.txt
# author Gabriel Rockefeller
# date   2011 June 13
# brief  Instructions for building device/test Makefile.
# note   © Copyright 2011 Los Alamos National Security, All rights reserved.
#------------------------------------------------------------------------------#
# $Id$

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

#set( CUDA_BUILD_CUBIN ON )
cuda_compile_ptx( generated_files ${cuda_sources} )

foreach( ptxfile ${generated_files} )
   get_filename_component( sh_ptxfile ${ptxfile} NAME_WE )
   string( REPLACE "cuda_compile_ptx_generated_" "" sh_ptxfile
      ${sh_ptxfile} )
   set( short_ptxfile "${sh_ptxfile}.ptx" )
   add_custom_target( RENAME_${sh_ptxfile} ALL
      COMMAND ${CMAKE_COMMAND} -E rename ${ptxfile} ${short_ptxfile}
      DEPENDS ${ptxfile} )
endforeach()

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
   Lib_dsxx
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


# set( test_deps
#    Lib_c4
#    Lib_dsxx
#    ${MPI_LIBRARIES}
#    )

# # Add tests
# add_parallel_tests(
#    SOURCES  "${test_sources}"
#    PE_LIST  "1;2;4"
#    DEPS     "${test_deps}"
#    MPIFLAGS "-q"
#    RESOURCE_LOCK "singleton"
#    )


