#-----------------------------*-cmake-*----------------------------------------#
# file   diagnostics/qt/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2015 June 26
# brief  Generate build project files for diagnostics/qt
# note   Copyright (C) 2016-2017, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
cmake_minimum_required( VERSION 3.9.0 )
project( diagnostics_qt CXX )

# ---------------------------------------------------------------------------- #
# Special options for Qt applications
# ---------------------------------------------------------------------------- #

# Find generated files in the corresponding build directories:
# 1. <file>.ui files generate ui_<file>.h headers in the build directory via the
#    qt5_wrap_ui macro.
# 2. <file>.qrc files generate qrc_<file>.cpp in the build directory via the
#    qt5_add_resources macro.

# Instruct CMake to run moc automatically when needed (only for subdirectories
# that need Qt)
set(CMAKE_AUTOMOC ON)

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB sources   *.cc *.cpp )
file( GLOB headers   *.hh *.h )
file( GLOB ui_files  *.ui )
#file( GLOB resources *.qrc )

# Make the header files available in the IDE.
if( MSVC_IDE OR ${CMAKE_GENERATOR} MATCHES Xcode )
   list( APPEND sources ${headers} )
endif()

# Use moc to convert *.ui files into ui_*.h files:
qt5_wrap_ui( ui_headers ${ui_files} )
list( APPEND sources ${ui_headers} )

# use rcc to convert *.qrc files into qrc_*.cpp files:
qt5_add_resources( qrc_sources ${resources} )
list( APPEND sources ${qrc_sources} )

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories(
   ${PROJECT_SOURCE_DIR}  # sources
   ${PROJECT_BINARY_DIR}  # qt generated .h files
   ${draco_src_dir_SOURCE_DIR} # ds++ and c4 headers
   ${dsxx_BINARY_DIR}     # ds++/config.h
   ${c4_BINARY_DIR}       # c4/config.h
   ${MPI_CXX_INCLUDE_PATH}
   )

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

# List Qt5::Widgets as a dependencies but not a VENDOR_LIB to prevent DBS from
# listing this library as a DRACO_TPL_LIBRARY.  It is not intended for export
# into other projects.
add_component_executable(
  TARGET      Exe_draco_info_gui
  TARGET_DEPS "Lib_diagnostics;Qt5::Widgets"
  SOURCES     "${sources}"
  PREFIX       Draco
  FOLDER       diagnostics
  NOCOMMANDWINDOW
  )

# Copy Qt dll files to build directory
# copy_dll_link_libraries_to_build_dir( Exe_draco_info_gui )
if (WIN32)
  file( GLOB qtdlls ${QTDIR}/bin/*.dll )
  add_custom_command( TARGET Exe_draco_info_gui POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${qtdlls}
            $<TARGET_FILE_DIR:Exe_draco_info_gui>
    COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:Lib_dsxx>
            $<TARGET_FILE_DIR:Exe_draco_info_gui>
    COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:Lib_c4>
            $<TARGET_FILE_DIR:Exe_draco_info_gui>
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:Lib_diagnostics>
            $<TARGET_FILE_DIR:Exe_draco_info_gui> )
endif()

# ---------------------------------------------------------------------------- #
# end diagnostics/qt/CMakeLists.txt
# ---------------------------------------------------------------------------- #
