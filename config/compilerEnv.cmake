#-----------------------------*-cmake-*----------------------------------------#
# file   config/compiler_env.cmake
# brief  Default CMake build parameters
# note   Copyright (C) 2010-2012 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

# Global

include( FeatureSummary )

# Library type to build
# Linux: STATIC is a lib<XXX>.a
#        SHARED is a lib<XXX>.so (requires rpath or .so found in $LD_LIBRARY_PATH
# MSVC : STATIC is <XXX>.lib
#        SHARED is <XXX>.dll (requires dll to be in $PATH or in same directory as exe).
if( NOT DEFINED DRACO_LIBRARY_TYPE )
   set( DRACO_LIBRARY_TYPE "SHARED" )
endif()
set( DRACO_LIBRARY_TYPE "${DRACO_LIBRARY_TYPE}" CACHE STRING 
	"Keyword for creating new libraries (STATIC or SHARED).")
# Provide a constrained drop down list in cmake-gui.
set_property( CACHE DRACO_LIBRARY_TYPE
   PROPERTY STRINGS SHARED STATIC)

if( EXISTS $ENV{PAPI_HOME} )
    set( HAVE_PAPI 1 CACHE BOOL "Is PAPI available on this machine?" )
    set( PAPI_INCLUDE $ENV{PAPI_INCLUDE} CACHE PATH 
       "PAPI headers at this location" )
    set( PAPI_LIBRARY $ENV{PAPI_LIBDIR}/libpapi.so CACHE FILEPATH
       "PAPI library." )
    if( NOT EXISTS ${PAPI_LIBRARY} )
       message( FATAL_ERROR "PAPI requested, but library not found.
    If on Turing, set PAPI_LIBDIR to correct path (module file is
    broken)." )
    endif()
    mark_as_advanced( PAPI_INCLUDE PAPI_LIBRARY )

    add_feature_info( HAVE_PAPI HAVE_PAPI 
       "Provide PAPI hardware counters if available." )
endif()


#----------------------------------------------------------------------#
# Macro to establish which runtime libraries to link against 
#
# Control link behavior for Run-Time Library.
# /MT - Causes your application to use the multithread, static
#       version of the run-time library. Defines _MT and causes
#       the compiler to place the library name LIBCMT.lib into the
#       .obj file so that the linker will use LIBCMT.lib to
#       resolve external symbols. 
# /MTd - Defines _DEBUG and _MT. This option also causes the
#       compiler to place the library name LIBCMTD.lib into the
#       .obj file so that the linker will use LIBCMTD.lib to
#       resolve external symbols. 
# /MD - Causes appliation to use the multithread and DLL specific
#       version of the run-time library.  Places MSVCRT.lib into
#       the .obj file.
#       Applications compiled with this option are statically
#       linked to MSVCRT.lib. This library provides a layer of
#       code that allows the linker to resolve external
#       references. The actual working code is contained in
#       MSVCR90.DLL, which must be available at run time to
#       applications linked with MSVCRT.lib. 
# /MD /D_STATIC_CPPLIB - applications link with the static
#       multithread Standard C++ Library (libcpmt.lib) instead of
#       the dynamic version (msvcprt.lib), but still links
#       dynamically to the main CRT via msvcrt.lib. 
# /MDd - Defines _DEBUG, _MT, and _DLL and causes your application
#       to use the debug multithread- and DLL-specific version of
#       the run-time library. It also causes the compiler to place
#       the library name MSVCRTD.lib into the .obj file. 
#----------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Setup compilers
#------------------------------------------------------------------------------#
macro(dbsSetupCompilers)

   # Bad platform
   if( NOT WIN32 AND NOT UNIX)
      message( FATAL_ERROR "Unsupported platform (not WIN32 and not UNIX )." )
   endif()  
   
   # shared or static libararies?
   if( ${DRACO_LIBRARY_TYPE} MATCHES "STATIC" )
      message(STATUS "Building static libraries.")
      set( MD_or_MT "MD" )
      set( DRACO_SHARED_LIBS 0 )
   elseif( ${DRACO_LIBRARY_TYPE} MATCHES "SHARED" )
      message(STATUS "Building shared libraries.")
      set( MD_or_MT "MD" )
      # This CPP symbol is used by config.h to signal if we are need to add 
      # declspec(dllimport) or declspec(dllexport) for MSVC.
      set( DRACO_SHARED_LIBS 1 )
      mark_as_advanced(DRACO_SHARED_LIBS)
   else()
      message( FATAL_ERROR "DRACO_LIBRARY_TYPE must be set to either STATIC or SHARED.")
   endif()
   set( DRACO_SHARED_LIBS ${DRACO_SHARED_LIBS} CACHE STRING 
      "This CPP symbol is used by config.h to signal if we are need to add declspec(dllimport) or declspec(dllexport) for MSVC." )
   
   ##---------------------------------------------------------------------------##
   ## Check for OpenMP
   
   set( USE_OPENMP OFF )
   include( CheckIncludeFiles )
   check_include_files( omp.h HAVE_OMP_H )
   if( ${HAVE_OMP_H} )
      set( USE_OPENMP ON )
   endif()   
   option( USE_OPENMP "Turn on OpenMP features?" ${USE_OPENMP} )

   ##---------------------------------------------------------------------------##

   set( gen_comp_env_set 1 )
endmacro()

#------------------------------------------------------------------------------#
# Setup C++ Compiler
#------------------------------------------------------------------------------#
macro(dbsSetupCxx)
   
   if( NOT gen_comp_env_set STREQUAL 1 )
      dbsSetupCompilers()
   endif()
   
   if( ${CMAKE_CXX_COMPILER} MATCHES "tau_cxx.sh" )
      execute_process(
         COMMAND ${CMAKE_CXX_COMPILER} -tau:showcompiler
         OUTPUT_VARIABLE my_cxx_compiler
         )
   else()
      set( my_cxx_compiler ${CMAKE_CXX_COMPILER} )
   endif()
   
   string( REGEX REPLACE ".*([0-9]).([0-9]).([0-9]).*" "\\1"
      DBS_CXX_COMPILER_VER_MAJOR "${CMAKE_CXX_COMPILER_VERSION}" )
   string( REGEX REPLACE ".*([0-9]).([0-9]).([0-9]).*" "\\2"
      DBS_CXX_COMPILER_VER_MINOR "${CMAKE_CXX_COMPILER_VERSION}" )
   
   if( ${my_cxx_compiler} MATCHES "cl" )
      include( windows-cl )
   elseif( ${my_cxx_compiler} MATCHES "ppu-g[+][+]" )
      include( unix-ppu )
   elseif( ${my_cxx_compiler} MATCHES "spu-g[+][+]" )
      include( unix-spu )
   elseif( ${my_cxx_compiler} MATCHES "icpc" )
      include( unix-intel )
   elseif( ${my_cxx_compiler} MATCHES "pgCC" )
      include( unix-pgi )
   elseif( ${my_cxx_compiler} MATCHES "xt-asyncpe" ) # Ceilo (catamount/pgi)
      include( unix-pgi )
   elseif( ${my_cxx_compiler} MATCHES "[cg][+]+" )
      include( unix-g++ )
   else( ${my_cxx_compiler} MATCHES "cl" )
      message( FATAL_ERROR "Build system does not support CXX=${my_cxx_compiler}" )
   endif( ${my_cxx_compiler} MATCHES "cl" )

endmacro()

#------------------------------------------------------------------------------#
# Setup Fortran Compiler
#
# Use:
#    include( compilerEnv )
#    dbsSetupFortran( [QUIET] )
#
# Returns:
#    BUILD_SHARED_LIBS - bool
#    CMAKE_Fortran_COMPILER - fullpath
#    CMAKE_Fortran_FLAGS
#    CMAKE_Fortran_FLAGS_DEBUG
#    CMAKE_Fortran_FLAGS_RELEASE
#    CMAKE_Fortran_FLAGS_RELWITHDEBINFO
#    CMAKE_Fortran_FLAGS_MINSIZEREL
#    ENABLE_SINGLE_PRECISION - bool
#    DBS_FLOAT_PRECISION     - string (config.h)
#    PRECISION_DOUBLE | PRECISION_SINGLE - bool
# 
#------------------------------------------------------------------------------#
macro(dbsSetupFortran)
   
   if( ${CMAKE_Fortran_COMPILER} MATCHES "gfortran" )
      include( unix-gfortran )
   elseif( ${CMAKE_Fortran_COMPILER} MATCHES "ifort" )
      include( unix-ifort )
   elseif( ${CMAKE_Fortran_COMPILER} MATCHES "pgf9[05]" )
      include( unix-pgf90 )
   elseif( ${my_cxx_compiler} MATCHES "xt-asyncpe" ) # Ceilo (catamount/pgi)
      include( unix-pgf90 )
   elseif( ${my_cxx_compiler} MATCHES "pgf90" )
      include( unix-pgf90 )
   else()
      message( FATAL_ERROR "Build system does not support F90=${CMAKE_Fortran_COMPILER}" )
   endif()
   
endmacro()

#------------------------------------------------------------------------------#
# End config/compiler_env.cmake
#------------------------------------------------------------------------------#
