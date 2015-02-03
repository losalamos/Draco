#-----------------------------*-cmake-*----------------------------------------#
# file   ds++/bin/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2013 Nov 6
# brief  Provide.
# note   Copyright (C) 2010-2015, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id: CMakeLists.txt 7273 2013-10-23 18:43:16Z kellyt $
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.0.0)
project( dsxx_bin CXX )

if( WIN32 )
# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories(
   ${PROJECT_SOURCE_DIR}       # sources
   ${draco_src_dir_SOURCE_DIR} # ds++ headers
   ${dsxx_BINARY_DIR}          # ds++/config.h
   )

# ---------------------------------------------------------------------------- #
# Build binaries
# ---------------------------------------------------------------------------- #

add_executable( Exe_testCrashHandler ${PROJECT_SOURCE_DIR}/testCrashHandler.cc )
target_link_libraries( Exe_testCrashHandler Lib_dsxx  )
set_target_properties( Exe_testCrashHandler
   PROPERTIES 
      OUTPUT_NAME testCrashHandler 
      VS_KEYWORD  testCrashHandler
      FOLDER      dsxx_test
)

# copy dll to bin directory.
# get_target_property( Lib_dsxx_loc Lib_dsxx LOCATION )
set( Lib_dsxx_loc $<TARGET_FILE:Lib_dsxx> )
# Also grab the file with debug info
string( REPLACE ".dll" ".pdb" pdb_file ${Lib_dsxx_loc} )
add_custom_command( TARGET Exe_testCrashHandler
   POST_BUILD
   COMMAND ${CMAKE_COMMAND} -E copy_if_different ${Lib_dsxx_loc}
           ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
   COMMAND ${CMAKE_COMMAND} -E copy_if_different ${pdb_file}
           ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
   )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Exe_testCrashHandler DESTINATION bin )

endif() # WIN32