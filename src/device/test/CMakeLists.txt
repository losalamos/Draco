#-----------------------------*-cmake-*----------------------------------------#
# file   device/test/CMakeLists.txt
# author Gabriel Rockefeller
# date   2011 June 13
# brief  Instructions for building device/test Makefile.
# note   © Copyright 2011 Los Alamos National Security, All rights reserved.
#------------------------------------------------------------------------------#
# $Id$

project( device_test CXX )

# ---------------------------------------------------------------------------- #
# Build one of DACS_Device or GPU_Device
# ---------------------------------------------------------------------------- #

if( ENABLE_HET_CODE )
   include( DACS_Device.cmake )
elseif( EXISTS "${CUDA_CUDA_LIBRARY}" )
   include( GPU_Device.cmake )
endif()

