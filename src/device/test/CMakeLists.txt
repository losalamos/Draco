#-----------------------------*-cmake-*----------------------------------------#
# file   device/test/CMakeLists.txt
# brief  Instructions for building device/test Makefile.
# note   Copyright(C) 2016-2019 Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
project( device_test CXX ${CUDA_DBS_STRING} )

# ---------------------------------------------------------------------------- #
# Build GPU_Device
# ---------------------------------------------------------------------------- #
if( WITH_CUDA )
  include( GPU_Device.cmake )
endif()

# ---------------------------------------------------------------------------- #
# End device/test/CMakeLists.txt
# ---------------------------------------------------------------------------- #
