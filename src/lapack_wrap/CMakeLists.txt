#-----------------------------*-cmake-*----------------------------------------#
# file   config/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for lapack_wrap.
# note   Copyright (C) 2010-2012, Los Alamos National Security
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 2.6)
project( lapack_wrap CXX )

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #
if( NOT lapack_FOUND OR APPLE )
   if( NOT lapack_FOUND )
      message( STATUS
      "NOTE: Skipping lapack_wrap component because LAPACK is not found." )
   else()
      message( STATUS
      "NOTE: Skipping lapack_wrap component because SYSTEM is APPLE" )
   endif()
else()

if( MSVC )

   # message("CMakeAddFortranSubdirectory")
   # include(CMakeAddFortranSubdirectory)
   # cmake_add_fortran_subdirectory( f90proj 
       # PROJECT f90proj
       # ARCHIVE_DIR ${CURRENT_BINARY_DIR}/f90proj 
       # RUNTIME_DIR ${CURRENT_BINARY_DIR}/f90proj 
       # NO_EXTERNAL_INSTALL
       # )
       
   # cmake_add_fortran_subdirectory( fortran
      # PROJECT DGESV
      # ARCHIVE_DIR ${CMake_ARCHIVE_OUTPUT_DIRECTORY}
      # RUNTIME_DIR ${CMake_RUNTIME_OUTPUT_DIRECTORY}
      # LIBRARIES refblas dgesv # target libraries created
      # LINK_LIBRARIES  # link interface libraries
      # LINK_LIBS dgesv refblas  # dgesv needs refblas to link
      # CMake_COMMAND_LINE
      # -DCMake_ARCHIVE_OUTPUT_DIRECTORY=${CMake_ARCHIVE_OUTPUT_DIRECTORY}
      # -DCMake_RUNTIME_OUTPUT_DIRECTORY=${CMake_RUNTIME_OUTPUT_DIRECTORY}
      # NO_EXTERNAL_INSTALL
   # )    
    
    # The autofortran library linking does not work unless fortran is enabled.
    # include(CheckLanguage)
    # check_language(Fortran)
    # if( CMake_Fortran_Compiler )
       # enable_language(Fortran)
    # else()
       # message(STATUS "No native Fortran support.")
    # endif()
    
    set( file_content "
#ifndef FC_HEADER_INCLUDED
#define FC_HEADER_INCLUDED

/* Mangling for Fortran global symbols without underscores. */
#define FC_GLOBAL(name,NAME) name##_

/* Mangling for Fortran global symbols with underscores. */
#define FC_GLOBAL_(name,NAME) name##_

/* Mangling for Fortran module symbols without underscores. */
#define FC_MODULE(mod_name,name, mod_NAME,NAME) mod_name##_mp_##name##_

/* Mangling for Fortran module symbols with underscores. */
#define FC_MODULE_(mod_name,name, mod_NAME,NAME) mod_name##_mp_##name##_

#endif
")
    if( NOT EXISTS ${PROJECT_BINARY_DIR}/lapack_wrap/config.h )
      file( WRITE ${PROJECT_BINARY_DIR}/lapack_wrap/config.h ${file_content} )
    endif()
    
else()

# cmake --help-module FortranCInterface
# This module automatically detects the API by which C and Fortran
# languages interact.
include(FortranCInterface)

# A function is provided to generate a C header file containing macros
# to mangle symbol names: 
#
# FortranCInterface_HEADER(<file>
#                          [MACRO_NAMESPACE <macro-ns>]
#                          [SYMBOL_NAMESPACE <ns>]
#                          [SYMBOLS [<module>:]<function> ...])
FortranCInterface_HEADER( lapack_wrap/config.h 
    MACRO_NAMESPACE "FC_"
    SYMBOLS scopy dcopy saxpy daxpy sscal dscal sdot ddot snrm2 dnrm2
    SYMBOL_NAMESPACE "FC_"
    )
# This creates a "lapack_wrap.h" header that defines mangling macros
# FC_GLOBAL(), FC_GLOBAL_() FC_MODULE(), and FC_MODULE().
# define saxpy FC_GLOBAL(saxpy, SAXPY)

endif()

# Omit this package is the mangling was not found by
# FortranCInterface.cmake 
if( EXISTS ${PROJECT_BINARY_DIR}/lapack_wrap/config.h )
# if( FortranCInterface_GLOBAL_FOUND AND FortranCInterface_MODULE_FOUND )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB sources *.cc )
file( GLOB headers *.hh )
list( APPEND headers ${PROJECT_BINARY_DIR}/lapack_wrap/config.h )

# Make the header files available in the IDE.
if( MSVC_IDE OR ${CMAKE_GENERATOR} MATCHES Xcode )
   list( APPEND sources ${headers} )
endif()

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories( 
   ${PROJECT_SOURCE_DIR}     # sources
   ${PROJECT_BINARY_DIR}     # lapack_wrap/config.h
   ${PROJECT_SOURCE_DIR}/..  # ds++ header files
   ${dsxx_BINARY_DIR}        # ds++/config.h
) 

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( FILES ${headers} DESTINATION include/lapack_wrap )

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

# --------------------------------------------------------------------------- #
# FortranCInterface_GLOBAL_FOUND AND FortranCInterface_MODULE_FOUND 
else()
   message("   lapack_wrap - skipped because FortranCInterface failed.")
endif()

# --------------------------------------------------------------------------- #
# lapack_FOUND?
else()
   message("   lapack_wrap - skipped because lapack was not found.")
endif() # NOT lapack_FOUND
