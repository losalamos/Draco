# $Id$

cmake_minimum_required(VERSION 2.6)
project( cdi_ipcress CXX )

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

#file( GLOB sources *.cc )
#file( GLOB explicit_instantiations *_pt.cc )
#file( GLOB headers *.hh )
#list( APPEND headers ${PROJECT_BINARY_DIR}/cdi_ipcress/config.h )

set( sources
   IpcressFile.cc
   IpcressDataTable.cc
   IpcressGrayOpacity.cc
   IpcressMultigroupOpacity.cc
)
set( headers
   IpcressFile.hh
   IpcressFile.t.hh
   IpcressDataTable.hh
   IpcressMultigroupOpacity.hh
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
#   ${PROJECT_BINARY_DIR}      # cdi_ipcress/config.h
   ${dsxx_BINARY_DIR}         # ds++/config.h
) 

# setup binary directory
#file( MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/bin )

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library( Lib_cdi_ipcress ${PROJECT_NAME} "${sources}" )
add_dependencies(      Lib_cdi_ipcress Lib_cdi Lib_dsxx )


#add_executable( Exe_Ipcress_Interpreter ${sources} )
#target_link_libraries( Exe_Ipcress_Interpreter 
#                       Lib_cdi_ipcress Lib_cdi Lib_dsxx 
#                       ${IPCRESS_LIBRARY} ${IPCRESS_GFORTRAN_LIBRARY})
#set_target_properties( Exe_Ipcress_Interpreter 
#                       PROPERTIES OUTPUT_NAME bin/Ipcress_Interpreter )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

#install( TARGETS Lib_cdi_ipcress DESTINATION lib )
#install( FILES ${headers} DESTINATION include/cdi_ipcress )
#install( TARGETS Exe_Ipcress_Interpreter DESTINATION bin )

# ---------------------------------------------------------------------------- #
# Unit tests
# ---------------------------------------------------------------------------- #

if( BUILD_TESTING )
   add_subdirectory( test )
endif()   

# ---------------------------------------------------------------------------- #
# Autodoc
# ---------------------------------------------------------------------------- #

#process_autodoc_pages()


