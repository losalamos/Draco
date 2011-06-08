# $Id$

cmake_minimum_required(VERSION 2.6)
project( lapack_wrap CXX )

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #
if( NOT LAPACK_FOUND )
   message( STATUS
      "NOTE: Skipping lapack_wrap component because LAPACK is not found." )
else()

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
# Build package library
# ---------------------------------------------------------------------------- #

# No sources!
# add_component_library( Lib_lapack_wrap ${PROJECT_NAME} "${sources}"  )

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

endif( NOT LAPACK_FOUND )

