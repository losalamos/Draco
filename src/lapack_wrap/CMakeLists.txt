# $Id$

cmake_minimum_required(VERSION 2.6)
project( lapack_wrap CXX )

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #
if( LAPACK_FOUND )
   # include_directories( ${lapack_INCLUDE_DIRS} )
   # include( CheckIncludeFiles )
   # usage: CHECK_INCLUDE_FILES (<header> <RESULT_VARIABLE> )
   # check_include_files( fblaswr.h HAVE_FBLASWR_H )
   # check_symbol_exists( scopy_ fblaswr.h HAVE_SCOPY )
   # if( HAVE_SCOPY )
      # set( SCOPY scopy_ )
   # endif()
else()
   find_package( LAPACK REQUIRED )
endif()
set( SCOPY scopy_ )
set( DCOPY dcopy_ )
set( SSCAL sscal_ )
set( SSCAL sscal_ )
set( DSCAL dscal_ )
set( SDOT  sdot_ )
set( DDOT  ddot_ )
set( DAXPY daxpy_ )
set( SAXPY saxpy_ )
set( SNRM2 snrm2_ )
set( DNRM2 dnrm2_ )
configure_file( config.h.in ${PROJECT_BINARY_DIR}/lapack_wrap/config.h )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB sources *.cc )
file( GLOB headers *.hh )
list( APPEND headers ${PROJECT_BINARY_DIR}/lapack_wrap/config.h )

# Make the header files available in the IDE.
if( MSVC_IDE )
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
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library( Lib_lapack_wrap ${PROJECT_NAME} "${sources}"  )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_lapack_wrap DESTINATION lib )
install( FILES ${headers} DESTINATION include/lapack_wrap )

# ---------------------------------------------------------------------------- #
# Unit tests
# ---------------------------------------------------------------------------- #

if( DRACO_BUILD_TESTS )
   add_subdirectory( test )
endif()   
