#-----------------------------*-cmake-*----------------------------------------#
# file   cdi_ipcress/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for cdi_ipcress.
# note   Copyright (C) 2010-2013, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#
cmake_minimum_required( VERSION 2.8 )
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

# Make the header files available in the IDE.
if( MSVC_IDE OR ${CMAKE_GENERATOR} MATCHES Xcode )
   list( APPEND sources ${headers} )
endif()

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories( 
   ${PROJECT_SOURCE_DIR}      # sources
   ${PROJECT_SOURCE_DIR}/..   # header files
   ${dsxx_BINARY_DIR}         # ds++/config.h
) 

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library( 
   TARGET       Lib_cdi_ipcress
   TARGET_DEPS  "Lib_cdi;Lib_c4"
   LIBRARY_NAME ${PROJECT_NAME}
   SOURCES      "${sources}" )

add_executable( Exe_Ipcress_Interpreter IpcressInterpreter.cc )
target_link_libraries( 
   Exe_Ipcress_Interpreter 
   Lib_cdi_ipcress Lib_cdi Lib_dsxx 
   ${IPCRESS_LIBRARY} ${IPCRESS_GFORTRAN_LIBRARY})
set_target_properties( 
   Exe_Ipcress_Interpreter 
   PROPERTIES OUTPUT_NAME Ipcress_Interpreter )
if( WIN32 )   
   # Add a post-build command to copy each dll into the test directory.
   foreach( lib Lib_cdi_ipcress;Lib_cdi;Lib_dsxx )
      unset( ${comp_target}_loc )
      get_target_property( ${comp_target}_loc ${lib} LOCATION )
      # Also grab the file with debug info
      string( REPLACE ".dll" ".pdb" pdb_file ${${comp_target}_loc} )
      add_custom_command( TARGET Exe_Ipcress_Interpreter 
         POST_BUILD
         COMMAND ${CMAKE_COMMAND} -E copy_if_different ${${comp_target}_loc} 
                 ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
         COMMAND ${CMAKE_COMMAND} -E copy_if_different ${pdb_file} 
                 ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
         )
   endforeach()
endif()
      
      

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #
set( CMAKE_INSTALL_DO_STRIP NO )
install( TARGETS Lib_cdi_ipcress EXPORT draco-targets DESTINATION lib )
install( FILES ${headers} DESTINATION include/cdi_ipcress )
install( TARGETS Exe_Ipcress_Interpreter DESTINATION bin )

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



