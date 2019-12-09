#-----------------------------*-cmake-*----------------------------------------#
# file   plot2D/test/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for plot2D/test.
# note   Copyright (C) 2016-2019, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
project( plot2D_test CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB test_sources *.cc )
file( GLOB bench_files bench/*.agr bench/*.par )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

set( test_deps Lib_plot2D ${Grace_LIBRARY} )

# This test requires numdiff
find_program( NUMDIFF
              numdiff
              HINTS ${VENDOR_DIR}
            )
if ( NOT ${NUMDIFF} MATCHES NOTFOUND )
   add_scalar_tests(
      SOURCES "${test_sources}"
      DEPS    "${test_deps}" )
   ##------------------------------------------------------------------------##
   ## This test doesn't behave correctly under valgrind (the pause doesn't work
   ## and that causes everything else to fail).
   ## ------------------------------------------------------------------------##
   set_tests_properties( plot2D_tstPlot2D
      PROPERTIES LABELS "nomemcheck" )
endif()

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

#------------------------------------------------------------------------------#
# End plot2D/test/CMakeLists.txt
#------------------------------------------------------------------------------#
