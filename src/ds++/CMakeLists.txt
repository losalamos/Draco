#-----------------------------*-cmake-*----------------------------------------#
# file   ds++/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 Sep 24
# brief  Generate build project files for ds++.
# note   Copyright (C) 2016-2019, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.9.0)
project( dsxx CXX )

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #

# Review config/platform_checks.cmake for many of these settings.

# Choose platforms for fpe_trap
if( ${CMAKE_SYSTEM} MATCHES "OSF" )
   set( FPETRAP_OSF_ALPHA 1 )
elseif( ${CMAKE_SYSTEM} MATCHES "Linux" )
   # If the OS is Linux, fpe_trap should work as long as feenableexcept is
   # available.  That is, it should work for both GNU and Intel compiler suites.
   set( FPETRAP_LINUX_X86 1 )

   # Turn fpe_trap off for PGI compilers.
   if( ${CMAKE_C_COMPILER} MATCHES pgcc )
      set( FPETRAP_LINUX_X86 0 )
   endif()
elseif( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
   if( ${_CMAKE_OSX_MACHINE} MATCHES "ppc" )
      set( FPETRAP_DARWIN_PPC 1 )
   else()
      set( FPETRAP_DARWIN_INTEL 1 )
   endif()

elseif( ${CMAKE_SYSTEM_NAME} MATCHES "Windows" )
   set( FPETRAP_WINDOWS_X86 1 )
endif()

# Choose settings for DRACO_DIAGNOSTICS
# Uses a bit mask
# 001 - basic print messages (negligible performance hit)
# 010 - detailed print messages (some performance hit or intrusive
#       output)
# 100 - additional computation required for diagnostic messages
#       (including fpe_trap)
set( DRACO_DIAGNOSTICS "0"
   CACHE STRING "Diagnostics verbosity level (0-7). Default is 0." )
# provide a constrained drop down menu for cmake-gui
set_property( CACHE DRACO_DIAGNOSTICS PROPERTY STRINGS 0 1 2 3 4 5 6 7 )

# Provide an enum from CMAKE_CXX_COMPILER_VERSION
if( NOT "${CMAKE_CXX_COMPILER_VERSION}x" STREQUAL "x" )
  string( REPLACE "." "" CMAKE_CXX_COMPILER_VERSION_ENUM
    ${CMAKE_CXX_COMPILER_VERSION} )
endif()

# Is this a KNL
if( CRAY_PE AND "$ENV{CRAY_CPU_TARGET}" STREQUAL "mic-knl" )
  set( draco_isKNL ON )
endif()

# Information about the current machine
cmake_host_system_information(RESULT COMPILED_HOSTNAME QUERY HOSTNAME)
cmake_host_system_information(RESULT COMPILED_FQDN QUERY FQDN)
cmake_host_system_information(RESULT COMPILED_OS_NAME QUERY OS_NAME)
cmake_host_system_information(RESULT COMPILED_OS_RELEASE QUERY OS_RELEASE)
cmake_host_system_information(RESULT COMPILED_OS_VERSION QUERY OS_VERSION)
cmake_host_system_information(RESULT COMPILED_OS_PLATFORM QUERY OS_PLATFORM)
cmake_host_system_information(RESULT draco_is64bit QUERY IS_64BIT)
cmake_host_system_information(RESULT HAVE_HARDWARE_FPU QUERY HAS_FPU)
cmake_host_system_information(RESULT HAVE_HARDWARE_MMX QUERY HAS_MMX)
cmake_host_system_information(RESULT HAVE_HARDWARE_MMX_PLUS QUERY HAS_MMX_PLUS)
cmake_host_system_information(RESULT HAVE_HARDWARE_SSE QUERY HAS_SSE)
cmake_host_system_information(RESULT HAVE_HARDWARE_SSE2 QUERY HAS_SSE2)
cmake_host_system_information(RESULT HAVE_HARDWARE_SSE_FP QUERY HAS_SSE_FP)
cmake_host_system_information(RESULT HAVE_HARDWARE_SSE_MMX QUERY HAS_SSE_MMX)

# Create ds++/config.h
configure_file( config.h.in ${PROJECT_BINARY_DIR}/ds++/config.h )

# Generate the dll_declspec.h file that defines DLL_PUBLIC_<pkg> CPP macros.
include(generate_dll_declspec)
generate_dll_declspec( "ds++" "c4;cdi;cdi_test;cdi_analytic;cdi_eospac;cdi_ipcress;cdi_ipcress_test;device;diagnostics;dsxx;fit;FC_Derived_Type;linear;mesh;meshReaders;mesh_element;min;norms;ode;parser;quadrature;quadrature_test;rng;roots;RTT_Format_Reader;special_functions;timestep;timestep_test;units;viz" )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB template_implementations *.t.hh *.i.hh )
file( GLOB sources *.cc )
file( GLOB explicit_instantiations *_pt.cc )
file( GLOB dsxx_headers *.hh )
list( APPEND dsxx_headers
   ${PROJECT_BINARY_DIR}/ds++/config.h
   ${PROJECT_BINARY_DIR}/ds++/dll_declspec.h )
set( dsxx_terminal_headers
   ${PROJECT_SOURCE_DIR}/terminal/terminal.h
   ${PROJECT_SOURCE_DIR}/terminal/terminal_base.h )
list( REMOVE_ITEM dsxx_headers ${template_implementations} )

if( ${CMAKE_GENERATOR} MATCHES "Xcode" )
  # Xcode does not support per file, per configuration compile definitions, so
  # we add this compile definition for all files
  add_compile_definitions( BUILD_TYPE="$<CONFIG>" )
else()
  set_source_files_properties( Release.cc PROPERTIES
    COMPILE_DEFINITIONS BUILD_TYPE="$<CONFIG>" )
endif()

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library(
   TARGET       Lib_dsxx
   LIBRARY_NAME ds++
   SOURCES      "${sources}"
   HEADERS      "${dsxx_headers};${dsxx_terminal_headers}" )
target_link_libraries( Lib_dsxx ${Lib_win_winsock} )
target_include_directories( Lib_dsxx
  # source directory or install location
  PUBLIC
    $<BUILD_INTERFACE:${Draco_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include>
  # generated include directive files (config.h)
  PUBLIC $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}> )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_dsxx EXPORT draco-targets DESTINATION ${DBSCFGDIR}lib )
install( FILES ${dsxx_headers} ${template_implementations} DESTINATION
  ${DBSCFGDIR}include/ds++ )
install( FILES
    ${dsxx_terminal_headers}
    ${PROJECT_SOURCE_DIR}/terminal/LICENSE
    ${PROJECT_SOURCE_DIR}/terminal/README.md
  DESTINATION ${DBSCFGDIR}include/terminal )

# ---------------------------------------------------------------------------- #
# Unit tests
# ---------------------------------------------------------------------------- #

if( BUILD_TESTING )
   add_subdirectory( test )
endif()
if( WIN32 )
  add_subdirectory( bin )
endif()

# ---------------------------------------------------------------------------- #
# Autodoc
# ---------------------------------------------------------------------------- #

process_autodoc_pages()
