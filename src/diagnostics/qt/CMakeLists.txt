#-----------------------------*-cmake-*----------------------------------------#
# file   diagnostics/qt/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2014 June 26
# brief  Generate build project files for diagnostics/qt
# note   Copyright (C) 2010-2014, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#
cmake_minimum_required( VERSION 3.0.0 )
project( diagnostics_qt CXX )

##---------------------------------------------------------------------------##
# Special options for Qt applications
##---------------------------------------------------------------------------##

# Find generated files in the corresponding build directories:
# 1. <file>.ui files generate ui_<file>.h headers in the build
#    directory via the qt5_wrap_ui macro.
# 2. <file>.qrc files generate qrc_<file>.cpp in the build directory
#    via the qt5_add_resources macro.

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

add_executable(        Exe_draco_info_gui ${sources} )
target_link_libraries( Exe_draco_info_gui 
  Lib_diagnostics
  Qt5::Widgets )
set_target_properties( Exe_draco_info_gui
  PROPERTIES OUTPUT_NAME draco_info-gui
  PROJECT_LABEL draco_info-gui
  FOLDER        diagnostics )

# add_component_library( 
#    TARGET       Lib_diagnostics
#    TARGET_DEPS  "Lib_c4"
#    LIBRARY_NAME ${PROJECT_NAME}
#    SOURCES      "${sources}" )

# add_executable( Exe_draco_info ${PROJECT_SOURCE_DIR}/draco_info.cc )
# target_link_libraries( Exe_draco_info Lib_dsxx  )
# set_target_properties( Exe_draco_info
#    PROPERTIES OUTPUT_NAME draco_info 
#    PROJECT_LABEL draco_info
#    FOLDER        diagnostics )

# Copy Qt dll files to build directory
if (WIN32)
  file( GLOB qtdlls ${QTDIR}/bin/*.dll )
  foreach( qtdll ${qtdlls} )
    add_custom_command( TARGET Exe_draco_info_gui POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${qtdll} 
                               $<TARGET_FILE_DIR:Exe_draco_info_gui> )
  endforeach()
  
  add_custom_command( TARGET Exe_draco_info_gui POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:Lib_dsxx> $<TARGET_FILE_DIR:Exe_draco_info_gui>
    COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:Lib_c4> $<TARGET_FILE_DIR:Exe_draco_info_gui> 
    COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:Lib_diagnostics> $<TARGET_FILE_DIR:Exe_draco_info_gui> )
endif()

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

# install( TARGETS Lib_diagnostics
#    EXPORT draco-targets
#    DESTINATION lib )

# # ??? Do I need to put draco_info into a bin/ directory for RPATH to
# # work correctly after install?
# install( TARGETS Exe_draco_info DESTINATION bin )
# install( FILES ${headers}  DESTINATION include/diagnostics )

# ---------------------------------------------------------------------------- #
# Unit tests
# ---------------------------------------------------------------------------- #

# if( BUILD_TESTING )
#    add_subdirectory( test )
# endif()   

# ---------------------------------------------------------------------------- #
# Autodoc
# ---------------------------------------------------------------------------- #

# process_autodoc_pages()
