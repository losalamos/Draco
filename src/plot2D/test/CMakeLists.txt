#=============================================================================
# Draco - Support libraries for transport solvers.
# Copyright 2010 LANS, LLC.
#=============================================================================
# $Id$
#=============================================================================

project( plot2D_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

#set( test_lib_sources 
#	${PROJECT_SOURCE_DIR}/test_utils.cc 
#	${PROJECT_SOURCE_DIR}/test_utils.hh)
file( GLOB test_sources *.cc )
file( GLOB bench_files bench/*.agr bench/*.par )
#list( REMOVE_ITEM test_sources ${test_lib_sources} )

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #
include_directories( 
   ${PROJECT_SOURCE_DIR}      # headers for tests
   ${PROJECT_SOURCE_DIR}/..   # headers for package
   ${PROJECT_BINARY_DIR}/..   # config.h
   )

# ---------------------------------------------------------------------------- #
# Build library for test directory
# ---------------------------------------------------------------------------- #

#add_component_library( Lib_plot2D_test ${PROJECT_NAME} "${test_lib_sources}" )
#add_dependencies( Lib_plot2D_test Lib_plot2D )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

set( test_deps 
   Lib_dsxx 
   ${Grace_LIBRARY} 
   )

add_scalar_tests(
   SOURCES "${test_sources}"
   DEPS    "${test_deps}" ) 

# ---------------------------------------------------------------------------- #
# Copy input files to test directory (binary_dir)
# ---------------------------------------------------------------------------- #

provide_aux_files( FILES "${bench_files}" )

if( WIN32 )
   provide_aux_files( FILES "${Grace_DLL_LIBRARIES}" )
endif()

##---------------------------------------------------------------------------##
## Extra 'make clean' files
##---------------------------------------------------------------------------##
set( extra_clean_files
   bench_plot1.part0
   bench_plot1.part1
   bench_plot2.part0
   bench_plot2.part1
   bench_plot3.part0
   bench_plot3.part1
   bench_plot4.part0
   bench_plot4.part1
   plot1.agr
   plot1.part0
   plot1.part1
   plot2.agr
   plot2.part0
   plot2.part1
   plot3.agr
   plot3.part0
   plot3.part1
   plot4.agr
   plot4.part0
   plot4.part1
   )
set_directory_properties(
   PROPERTIES
   ADDITIONAL_MAKE_CLEAN_FILES "${extra_clean_files}" )

##---------------------------------------------------------------------------##
## This test doesn't behave correctly under valgrind (the pause
## doesn't work and that causes everything else to fail).
##---------------------------------------------------------------------------##

set_tests_properties( plot2D_tstPlot2D
   PROPERTIES LABELS "nomemcheck" )
