#-----------------------------*-cmake-*----------------------------------------#
# file   FortranChecks/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for FortranChecks.
# note   Copyright (C) 2010-2014, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#
cmake_minimum_required( VERSION 2.8 )
project( FortranChecks CXX )

# When using CXX/Fortran cross language linking and function calling
# the CMake module FortranCInterface.cmake can be used to deal with
# name mangling issues.  However, Draco has chosen to use
# ISO_C_BINDING instead and this makes manual name demangling
# unnecessary.

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #
file( GLOB Csources *.cc )

# ---------------------------------------------------------------------------- #
# Build package libraries
# ---------------------------------------------------------------------------- #

add_component_library( 
   TARGET       Lib_FC_Derived_Type 
   LIBRARY_NAME FC_Derived_Type
   SOURCES      "${Csources}" 
   NOEXPORT )
   
set_target_properties( Lib_FC_Derived_Type 
   PROPERTIES FOLDER FortranChecks )

if( _LANGUAGES_ MATCHES Fortran OR MSVC )
   # For Unix, this will simply run add_subdirectory(f90sub).  For MVSE, 
   # an ExternalProject will be created that uses MinGW's gfortran via
   # Makefile to generate the library external to the MVSE GUI.

   # Try to use use gfortran with MVSE via add_fortran_subdirectory.
   # \todo Should this be moved to windows-cl.cmake?)
   if( DRACO_SHARED_LIBS )
      set( CMAKE_GNUtoMS ON )
   endif()
   # http://cmake.org/cmake/help/v2.8.11/cmake.html#module:CMakeAddFortranSubdirectory
   include(CMakeAddFortranSubdirectory)
   # cmake will put some run scripts in $build_dir/f90sub
   # Another directory for the build is generated at $build_dir/FortranChecks_f90sub_build-prefix
   cmake_add_fortran_subdirectory( f90sub  # directory name
      PROJECT FortranChecks_f90sub # project name in toplevel CMakeLists.txt
                                   # Creates target named 'FortranChecks_f90sub_build' 
      ARCHIVE_DIR f90sub      # .lib location relative to root binary tree
      RUNTIME_DIR f90sub      # .dll location relative to root binary tree
      LIBRARIES rtt_FC_f90sub # Sub project will create a library with this name.
                              # --> librtt_FC_f90.dll.
      NO_EXTERNAL_INSTALL
   )
   # The following files are generated:
   # $build_dir/f90sub/
   #                   librtt_FC_f90sub.dll
   #                   librtt_FC_f90sub.lib
   # target name is FC_f90sub
   if( MSVC )
      set_target_properties( FortranChecks_f90sub_build 
         PROPERTIES FOLDER FortranChecks )
   endif()    
else()
   message( STATUS "Skipping FortranChecks because Fortran not found in _LANGUAGES_." )
endif()

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #
# No need to install this library as it is not used by anyone else.  This 
# component is only used to test Fortran/C binding.

# ---------------------------------------------------------------------------- #
# Unit tests
# ---------------------------------------------------------------------------- #
if( BUILD_TESTING )
   add_subdirectory( test )
   if( _LANGUAGES_ MATCHES Fortran )
      add_subdirectory( ftest )
   endif()
   
   # The cmake_add_fortran_subdirectory does not currently support generating 
   # executables via MinGW's gfortran from within MSVC.
   # if( MSVC )   
      # unset( cmake_extra_commands )
      # foreach( dep Lib_FC_Derived_Type rtt_FC_f90sub )
        # set( ${dep}_loc $<TARGET_FILE:${dep}> )
        # list( APPEND cmake_extra_commands "-D${dep}=${${dep}_loc}" )
      # endforeach()
      # cmake_add_fortran_subdirectory( ftest
         # PROJECT FortranChecks_ftest 
         # ARCHIVE_DIR ftest 
         # RUNTIME_DIR ftest 
         # LIBRARIES Ut_FortranChecks_ftest_f90main_exe 
         # LINK_LIBRARIES 
            # LINK_LIBS Ut_FortranChecks_ftest_f90main_exe Lib_FC_Derived_Type rtt_FC_f90sub # Lib_dsxx
         # CMAKE_COMMAND_LINE "${cmake_extra_commands}"
         # NO_EXTERNAL_INSTALL
      # )
   # endif()
endif()   

# ---------------------------------------------------------------------------- #
# Autodoc
# ---------------------------------------------------------------------------- #
process_autodoc_pages()
