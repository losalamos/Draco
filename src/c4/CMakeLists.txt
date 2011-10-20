# $Id$

cmake_minimum_required(VERSION 2.6)
project( c4 CXX )

# ---------------------------------------------------------------------------- #
# Generate config.h (only occurs when cmake is run)
# ---------------------------------------------------------------------------- #

if( "${DRACO_C4}" MATCHES "MPI" )
  set( C4_MPI 1 CACHE STRING "Are we building an MPI aware version? (c4/config.h)"  )
  mark_as_advanced( C4_MPI )
else()
  set( C4_SCALAR 1 CACHE STRING "Are we building a scalar-only version (no mpi in c4/config.h)?" )
  mark_as_advanced( C4_MPI )
endif()
if( UNIX )
   set( C4_MPICMD "mpirun ${MPIEXEC_POSTFLAGS_STRING} -np " )
elseif( WIN32 )
   set( C4_MPICMD "mpiexec -np " )
elseif( OSF1 )
   set( C4_MPICMD "prun -n " )
elseif( APPLE )
   set( C4_MPICMD "mpirun -np " )
else() # AIX
   set( C4_MPICMD "poe -procs " )
endif()
if( ${CMAKE_SYSTEM_NAME} MATCHES "Catamount" )
   set( C4_MPICMD "aprun -n " )
endif()

if( EXISTS $ENV{PAPI_HOME} )
    set( HAVE_PAPI 1 )
    set( PAPI_INCLUDE $ENV{PAPI_INCLUDE} )
    set( PAPI_LIBRARY $ENV{PAPI_LIBDIR}/libpapi.so )
endif()

configure_file( config.h.in ${PROJECT_BINARY_DIR}/c4/config.h )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB template_implementations *.t.hh *.i.hh )
file( GLOB sources *.cc )
file( GLOB explicit_instantiations *_pt.cc )
file( GLOB headers *.hh *.h)
list( APPEND headers ${PROJECT_BINARY_DIR}/c4/config.h )
list( REMOVE_ITEM headers ${template_implementations} )

# Make the header files available in the IDE.
if( MSVC_IDE OR ${CMAKE_GENERATOR} MATCHES Xcode )
   list( APPEND sources ${headers} )
endif()

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories( ${PROJECT_SOURCE_DIR}       # sources
                     ${PROJECT_BINARY_DIR}       # c4/config.h
                     ${draco_src_dir_SOURCE_DIR} # ds++ header files
                     ${dsxx_BINARY_DIR}          # ds++/config.h
                     ${MPI_INCLUDE_PATH}
                     )

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library( Lib_c4 c4 "${sources}" ${PAPI_LIBRARY} )
add_dependencies( Lib_c4 Lib_dsxx )

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_c4 DESTINATION lib )
install( FILES ${headers} DESTINATION include/c4 )

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


