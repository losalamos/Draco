#=============================================================================
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
#=============================================================================
# $Id$
#=============================================================================

project( special_functions_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB test_sources *.cc )

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories( 
   ${PROJECT_SOURCE_DIR}      # headers for tests
   ${PROJECT_SOURCE_DIR}/..   # headers for package
   ${ode_BINARY_DIR}          # ode/config.h
   ${GSL_INCLUDE_DIRS}
   )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

# There is no Lib_special_functions_test so can't use the
# add_scalar_tests() macro. 

add_scalar_tests( SOURCES "${test_sources}" #)
   DEPS    "Lib_units;Lib_ode;Lib_dsxx;${GSL_LIBRARIES}" )

# ---------------------------------------------------------------------------- #
# Copy vendor dll files to test directory (binary_dir)
# ---------------------------------------------------------------------------- #

if( WIN32 )
   provide_aux_files( FILES "${GSL_DLL_LIBRARIES}" )
endif()

