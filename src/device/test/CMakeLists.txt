#-----------------------------*-cmake-*----------------------------------------#
# file   device/test/CMakeLists.txt
# brief  Instructions for building device/test Makefile.
# note   Copyright(C) 2011-2014 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$

project( device_test CXX )

# ---------------------------------------------------------------------------- #
# Build one of DACS_Device or GPU_Device
# ---------------------------------------------------------------------------- #

if( ENABLE_HET_CODE )
   include( DACS_Device.cmake )
elseif( USE_CUDA AND EXISTS "${CUDA_NVCC_EXECUTABLE}" )
   include( GPU_Device.cmake )
endif()

