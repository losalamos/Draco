#-----------------------------*-cmake-*----------------------------------------#
# file   device/test/CMakeLists.txt
# author Gabriel Rockefeller
# date   2011 June 13
# brief  Instructions for building device/test Makefile.
# note   © Copyright 2011 Los Alamos National Security, All rights reserved.
#------------------------------------------------------------------------------#
# $Id$

project( device_test CXX )

# If this is an x86 build, build the unit tests.  If this is a ppe
# build, build the ppe helper binary.  The ppe binary must be
# installed (i.e., available in PPE_BINDIR) for the x86 tests to pass.
#
# Another possible approach might have been to build this package only
# during x86 builds and then somehow just build the ppe helper binary
# using the ppe compiler (which must be available, since this is a
# het-only package).  In that scenario, the ppe binary wouldn't need
# to be installed; it would just be built in the test directory, next
# to the unit tests that use it.  CMake really prefers the "one build
# tree-one toolchain" model, though, and it's probably better not to
# fight the system.

if( NOT "${CMAKE_CXX_COMPILER}" MATCHES "[sp]pu-g[+][+]" )  # if x86 build

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

   file( GLOB test_sources *.cc )
   list( REMOVE_ITEM test_sources ${PROJECT_SOURCE_DIR}/dacs_device_ppe.cc )

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

   include_directories( 
      ${PROJECT_SOURCE_DIR}        # headers for tests
      ${PROJECT_SOURCE_DIR}/..     # headers for package
      ${PROJECT_BINARY_DIR}/..     # config.h
      )

# ---------------------------------------------------------------------------- #
# Build Unit tests
# ---------------------------------------------------------------------------- #

   set( test_deps
      Lib_c4
      Lib_dsxx
      ${MPI_LIBRARIES}
      )

   # Add tests
   add_parallel_tests(
      SOURCES  "${test_sources}"
      PE_LIST  "1;2;4"
      DEPS     "${test_deps}"
      MPIFLAGS "-q"
      )

else() # if ppe build

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

   include_directories(
      ${PROJECT_SOURCE_DIR}        # headers for tests
      ${PROJECT_SOURCE_DIR}/..     # headers for package
      ${PROJECT_BINARY_DIR}/..     # config.h
      ${draco_src_dir_SOURCE_DIR}  # ds++ header files
      ${dsxx_BINARY_DIR}           # ds++/config.h
      )

# ---------------------------------------------------------------------------- #
# Build helper binary
# ---------------------------------------------------------------------------- #

   add_executable( Exe_dacs_device_ppe dacs_device_ppe.cc )
   set_target_properties( Exe_dacs_device_ppe
      PROPERTIES
      OUTPUT_NAME dacs_device_ppe_exe
      LINK_FLAGS  "-Wl,-m,elf64ppc"
      )

   target_link_libraries( Exe_dacs_device_ppe
      Lib_dsxx
      /usr/lib64/libdacs_hybrid.so
      )

   install( TARGETS Exe_dacs_device_ppe DESTINATION bin )

endif()
