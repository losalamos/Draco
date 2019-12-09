#-----------------------------*-cmake-*----------------------------------------#
# file   cdi_ipcress/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for cdi_ipcress.
# note   Copyright (C) 2016-2019, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.9.0)
project( cdi_ipcress CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

set( sources
   IpcressFile.cc
   IpcressDataTable.cc
   IpcressGrayOpacity.cc
   IpcressMultigroupOpacity.cc
   IpcressOdfmgOpacity.cc
)
set( headers
   IpcressFile.hh
   IpcressMaterial.hh
   IpcressFile.t.hh
   IpcressDataTable.hh
   IpcressGrayOpacity.hh
   IpcressMultigroupOpacity.hh
   IpcressOdfmgOpacity.hh
)

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library(
  TARGET       Lib_cdi_ipcress
  TARGET_DEPS  "Lib_cdi"
  LIBRARY_NAME ${PROJECT_NAME}
  HEADERS      "${headers}"
  SOURCES      "${sources}" )


add_component_executable(
  TARGET      Exe_Ipcress_Interpreter
  TARGET_DEPS "Lib_cdi_ipcress"
  VENDOR_LIBS "${IPCRESS_LIBRARY};${IPCRESS_GFORTRAN_LIBRARY}"
  SOURCES     ${PROJECT_SOURCE_DIR}/IpcressInterpreter.cc
  PREFIX       Draco
  )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #
set( CMAKE_INSTALL_DO_STRIP NO )
install( TARGETS Lib_cdi_ipcress EXPORT draco-targets
  DESTINATION ${DBSCFGDIR}lib )
install( FILES ${headers} DESTINATION ${DBSCFGDIR}include/cdi_ipcress )
install( TARGETS Exe_Ipcress_Interpreter EXPORT draco-targets DESTINATION
  ${DBSCFGDIR}bin )

# ---------------------------------------------------------------------------- #
# Unit tests
# ---------------------------------------------------------------------------- #

if( BUILD_TESTING )
   add_subdirectory( test )
endif()

# ---------------------------------------------------------------------------- #
# Autodoc
# ---------------------------------------------------------------------------- #

process_autodoc_pages()

#------------------------------------------------------------------------------#
# End cdi_ipcress/CMakeLists.txt
#------------------------------------------------------------------------------#
