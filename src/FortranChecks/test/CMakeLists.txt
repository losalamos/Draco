# ---------------------------------------------------------------------------- #
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
# ---------------------------------------------------------------------------- #
# $Id: CMakeLists.txt 6136 2011-09-14 00:16:01Z kellyt $
# ---------------------------------------------------------------------------- #

project( FortranChecks_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB test_sources *.cc )

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories( 
   ${PROJECT_SOURCE_DIR}      # headers for tests
   ${dsxx_SOURCE_DIR}/..      # headers for ds++
   ${dsxx_BINARY_DIR}   )     # ds++/config.h

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

add_executable( cppmain ${test_sources} )
target_link_libraries( cppmain Lib_dsxx Lib_FortranChecks )


add_scalar_tests( 
   SOURCES "${test_sources}" 
   DEPS    "Lib_dsxx;Lib_FortranChecks" )

