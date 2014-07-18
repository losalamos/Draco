#-----------------------------*-cmake-*----------------------------------------#
# file   config/compiler_env.cmake
# brief  Default CMake build parameters
# note   Copyright (C) 2010-2013 Los Alamos National Security, LLC.
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
set_property( CACHE DRACO_LIBRARY_TYPE PROPERTY STRINGS SHARED STATIC)

# OpenMP default setup
if( NOT DEFINED USE_OPENMP )
   option( USE_OPENMP "Turn on OpenMP features?" ON )
endif()

# ----------------------------------------
# PAPI
# ----------------------------------------
if( EXISTS $ENV{PAPI_HOME} )
    set( HAVE_PAPI 1 CACHE BOOL "Is PAPI available on this machine?" )
    set( PAPI_INCLUDE $ENV{PAPI_INCLUDE} CACHE PATH 
       "PAPI headers at this location" )
    set( PAPI_LIBRARY $ENV{PAPI_LIBDIR}/libpapi.so CACHE FILEPATH
       "PAPI library." )
endif()
# PAPI 4.2 on CT uses a different setup.
if( $ENV{PAPI_VERSION} MATCHES "[45].[0-9].[0-9]")
    set( HAVE_PAPI 1 CACHE BOOL "Is PAPI available on this machine?" )
    string( REGEX REPLACE ".*[ ][-]I(.*)$" "\\1" PAPI_INCLUDE
       $ENV{PAPI_INCLUDE_OPTS} )
    string( REGEX REPLACE ".*[ ][-]L(.*)[ ].*" "\\1" PAPI_LIBDIR
       $ENV{PAPI_POST_LINK_OPTS} )
endif()
if( HAVE_PAPI )
    set( PAPI_INCLUDE ${PAPI_INCLUDE} CACHE PATH 
       "PAPI headers at this location" )
    set( PAPI_LIBRARY ${PAPI_LIBDIR}/libpapi.so CACHE FILEPATH
       "PAPI library." )
    if( NOT EXISTS ${PAPI_LIBRARY} )
       message( FATAL_ERROR "PAPI requested, but library not found.  Set PAPI_LIBDIR to correct path." )
    endif()
    mark_as_advanced( PAPI_INCLUDE PAPI_LIBRARY )
    add_feature_info( HAVE_PAPI HAVE_PAPI 
       "Provide PAPI hardware counters if available." )
endif()

# ------------------------------------------------------------------------------
# Identify machine and save name in ds++/config.h
# ------------------------------------------------------------------------------
site_name( SITENAME )
string( REGEX REPLACE "([A-z0-9]+).*" "\\1" SITENAME ${SITENAME} )
if( ${SITENAME} MATCHES "c[it]" )
  set( SITENAME "Cielito" )
elseif( ${SITENAME} MATCHES "ml[0-9]+" OR ${SITENAME} MATCHES "ml-fey" OR
        ${SITENAME} MATCHES "lu[0-9]+" OR ${SITENAME} MATCHES "lu-fey" )
  set( SITENAME "Moonlight" )
elseif( ${SITENAME} MATCHES "ccscs[0-9]+" )
  # do nothing (keep the fullname)
endif()
set( SITENAME ${SITENAME} CACHE "STRING" "Name of the current machine" FORCE)

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
   
   # Defaults for 1st pass:

   # shared or static libararies?
   if( ${DRACO_LIBRARY_TYPE} MATCHES "STATIC" )
      # message(STATUS "Building static libraries.")
      set( MD_or_MT "MD" )
      set( DRACO_SHARED_LIBS 0 )
      set( DRACO_LIBEXT ".a" )
   elseif( ${DRACO_LIBRARY_TYPE} MATCHES "SHARED" )
      # message(STATUS "Building shared libraries.")
      set( MD_or_MT "MD" )
      # This CPP symbol is used by config.h to signal if we are need to add 
      # declspec(dllimport) or declspec(dllexport) for MSVC.
      set( DRACO_SHARED_LIBS 1 )
      mark_as_advanced(DRACO_SHARED_LIBS)
      set( DRACO_LIBEXT ".so" )
   else()
      message( FATAL_ERROR "DRACO_LIBRARY_TYPE must be set to either STATIC or SHARED.")
   endif()
   set( DRACO_SHARED_LIBS ${DRACO_SHARED_LIBS} CACHE STRING 
      "This CPP symbol is used by config.h to signal if we are need to add declspec(dllimport) or declspec(dllexport) for MSVC." )
   
   ##---------------------------------------------------------------------------##
   ## Check for OpenMP

   #if( NOT gen_comp_env_set STREQUAL 1 )
   # set( USE_OPENMP OFF )
   if( USE_OPENMP )
      include( CheckIncludeFiles )
      check_include_files( omp.h HAVE_OMP_H )
      if( NOT ${HAVE_OMP_H} )
         set( USE_OPENMP OFF )
      endif()   
   endif()
   set( USE_OPENMP ${USE_OPENMP} CACHE BOOL "Turn on OpenMP features?" )
   # endif()

   ##---------------------------------------------------------------------------##

   # set( gen_comp_env_set 1 )

endmacro()

#------------------------------------------------------------------------------#
# Setup C++ Compiler
#------------------------------------------------------------------------------#
macro(dbsSetupCxx)
   
   # if( NOT gen_comp_env_set STREQUAL 1 )
      dbsSetupCompilers()
   # endif()
   
   # Deal with compiler wrappers
   if( ${CMAKE_CXX_COMPILER} MATCHES "tau_cxx.sh" )
      # When using the TAU profiling tool, the actual compiler vendor
      # is hidden under the tau_cxx.sh script.  Use the following
      # command to determine the actual compiler flavor before setting
      # compiler flags (end of this macro).
      execute_process(
         COMMAND ${CMAKE_CXX_COMPILER} -tau:showcompiler
         OUTPUT_VARIABLE my_cxx_compiler )
   elseif( ${CMAKE_CXX_COMPILER} MATCHES "xt-asyncpe" ) 
      # Ceilo (catamount) uses a wrapper script
      # /opt/cray/xt-asyncpe/5.06/bin/CC that masks the actual
      # compiler.  Use the following command to determine the actual
      # compiler flavor before setting compiler flags (end of this
      # macro).
      execute_process(
         COMMAND ${CMAKE_CXX_COMPILER} --version
         OUTPUT_VARIABLE my_cxx_compiler
         ERROR_QUIET )
      string( REGEX REPLACE "^(.*).Copyright.*" "\\1" 
         my_cxx_compiler ${my_cxx_compiler})
      # If a wrapper script is used, CMake will not have found the
      # compiler version...
      # icpc (ICC) 12.1.2 20111128
      # pgCC 11.10-0 64-bit target 
      # g++ (GCC) 4.6.2 20111026 (Cray Inc.)
      if( "x${CMAKE_CXX_COMPILER_VERSION}" STREQUAL "x" )
         string( REGEX REPLACE ".* ([0-9]+[.][0-9]+[.-][0-9]+).*" "\\1"
            CMAKE_CXX_COMPILER_VERSION ${my_cxx_compiler} )
      endif()
   else()
      set( my_cxx_compiler ${CMAKE_CXX_COMPILER} )
   endif()
   
   string( REGEX REPLACE "[^0-9]*([0-9]+).([0-9]+).([0-9]+).*" "\\1"
      DBS_CXX_COMPILER_VER_MAJOR "${CMAKE_CXX_COMPILER_VERSION}" )
   string( REGEX REPLACE "[^0-9]*([0-9]+).([0-9]+).([0-9]+).*" "\\2"
      DBS_CXX_COMPILER_VER_MINOR "${CMAKE_CXX_COMPILER_VERSION}" )

   set( DRACO_ENABLE_STRICT_ANSI OFF CACHE BOOL "use strict ANSI flags, C98" )

   # C99 support:
   option( DRACO_ENABLE_C99 "Support C99 features." ON )

   # C++11 support:
   option( DRACO_ENABLE_CXX11 "Support C++11 features." ON )
   if( ${my_cxx_compiler} MATCHES "clang" OR 
         ${my_cxx_compiler} MATCHES "llvm")
      include( apple-clang )
   elseif( ${my_cxx_compiler} MATCHES "cl" )
      include( windows-cl )
   elseif( ${my_cxx_compiler} MATCHES "icpc" )
      include( unix-intel )
   elseif( ${my_cxx_compiler} MATCHES "pgCC" )
      include( unix-pgi )
   elseif( ${my_cxx_compiler} MATCHES "xl" )
      set( DRACO_ENABLE_CXX11 ON )
      include( unix-xl )
   elseif( ${my_cxx_compiler} MATCHES "[cg][+]+" )
      include( unix-g++ )
   else( ${my_cxx_compiler} MATCHES "cl" )
      message( FATAL_ERROR "Build system does not support CXX=${my_cxx_compiler}" )
   endif()

   # Force possibly new values to cache
   set( DRACO_ENABLE_CXX11 ${DRACO_ENABLE_CXX11} CACHE BOOL 
      "Support C++11 features." FORCE )
   set( DRACO_ENABLE_C99 ${DRACO_ENABLE_C99} CACHE BOOL 
      "Support C99 features." FORCE )
   set( DRACO_ENABLE_STRICT_ANSI ${DRACO_ENABLE_STRICT_ANSI} CACHE
      BOOL "use strict ANSI flags, C98" FORCE )

   # Sanity check
   if( (DRACO_ENABLE_CXX11 AND DRACO_ENABLE_STRICT_ANSI) OR
         (DRACO_ENABLE_C99 AND DRACO_ENABLE_STRICT_ANSI) )
      message( FATAL_ERROR "Both DRACO_ENABLE_CXX11=ON and "
         "DRACO_ENABLE_STRICT_ANSI=ON, or DRACO_ENABLE_C99=ON and "
         "DRACO_ENABLE_STRICT_ANSI=ON, are not supported because "
         "STRICT_ANSI implied C++98 standard only." )
   endif()

   # To the greatest extent possible, installed versions of packages
   # should record the configuration options that were used when they
   # were built.  For preprocessor macros, this is usually
   # accomplished via #define directives in config.h files.  A
   # package's installed config.h file serves as both a record of
   # configuration options and a central location for macro
   # definitions that control features in the package.  Defining
   # macros via the -D command-line option to the preprocessor leaves
   # no record of configuration choices (except in a build log, which
   # may not be preserved with the installation).
   #
   # Unfortunately, there are cases where a particular macro must be
   # defined before some particular system header file is included, or
   # before any system header files are included.  In these
   # situations, using the config.h mechanism introduces sensitivity
   # to the order of header files, which can lead to brittleness;
   # defining project-wide language- or system-feature macros via -D,
   # using CMake's add_definitions command, is an acceptable
   # alternative.  Such definitions appear below.

   # Enable the definition of UINT64_C in stdint.h (required by
   # Random123).
   add_definitions(-D__STDC_CONSTANT_MACROS)
   set( CMAKE_REQUIRED_DEFINITIONS
      "${CMAKE_REQUIRED_DEFINITIONS} -D__STDC_CONSTANT_MACROS" )

   # Define _POSIX_C_SOURCE=200112 and _XOPEN_SOURCE=600, to enable
   # definitions conforming to POSIX.1-2001, POSIX.2, XPG4, SUSv2,
   # SUSv3, and C99.  See the feature_test_macros(7) man page for more
   # information.
   if( DRACO_ENABLE_C99 )
      add_definitions(-D_POSIX_C_SOURCE=200112 -D_XOPEN_SOURCE=600)
      set( CMAKE_REQUIRED_DEFINITIONS
         "${CMAKE_REQUIRED_DEFINITIONS} -D_POSIX_C_SOURCE=200112" )
      set( CMAKE_REQUIRED_DEFINITIONS
         "${CMAKE_REQUIRED_DEFINITIONS} -D_XOPEN_SOURCE=600")
      if ( APPLE ) 
          # Defining the above requires adding POSIX extensions,
          # otherwise, include ordering still goes wrong on Darwin,
          # (i.e., putting fstream before iostream causes problems)
          # see https://code.google.com/p/wmii/issues/detail?id=89
          add_definitions(-D_DARWIN_C_SOURCE)
          set( CMAKE_REQUIRED_DEFINITIONS
            "${CMAKE_REQUIRED_DEFINITIONS} -D_DARWIN_C_SOURCE ")
      endif() 
  endif()

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

   #if( NOT gen_comp_env_set STREQUAL 1 )
      dbsSetupCompilers()
   #endif()

   # Deal with comiler wrappers
   if( ${CMAKE_Fortran_COMPILER} MATCHES "xt-asyncpe" ) 
      # Ceilo (catamount) uses a wrapper script
      # /opt/cray/xt-asyncpe/5.06/bin/CC that masks the actual
      # compiler.  Use the following command to determine the actual
      # compiler flavor before setting compiler flags (end of this
      # macro).
      execute_process(
         COMMAND ${CMAKE_Fortran_COMPILER} --version
         # COMMAND ${CMAKE_Fortran_COMPILER} -V
         OUTPUT_VARIABLE my_fc_compiler
         ERROR_QUIET )
      string( REGEX REPLACE "^(.*).Copyright.*" "\\1" 
         my_fc_compiler ${my_fc_compiler})
   else()
      set( my_fc_compiler ${CMAKE_Fortran_COMPILER} )
   endif()

   if( ${my_fc_compiler} MATCHES "pgf9[05]" OR
         ${my_fc_compiler} MATCHES "pgfortran" )
      include( unix-pgf90 )
   elseif( ${my_fc_compiler} MATCHES "ifort" )
      include( unix-ifort )
   elseif( ${my_fc_compiler} MATCHES "xl" )
      include( unix-xlf )
   elseif( ${my_fc_compiler} MATCHES "gfortran" )
      include( unix-gfortran )
   else()
      message( FATAL_ERROR "Build system does not support F90=${my_fc_compiler}" )
   endif()

endmacro()

##---------------------------------------------------------------------------##
## Toggle a compiler flag based on a bool
##
## Examples:
##   toggle_compiler_flag( USE_OPENMP         "-fopenmp"   "C;CXX;EXE_LINKER" "")
##   toggle_compiler_flag( DRACO_ENABLE_CXX11 "-std=c++0x" "CXX" "DEBUG")
##---------------------------------------------------------------------------##
macro( toggle_compiler_flag switch compiler_flag
      compiler_flag_var_names build_modes )

   # generate names that are safe for CMake RegEx MATCHES commands
   string(REPLACE "+" "x" safe_compiler_flag ${compiler_flag})      

   # Loop over types of variables to check: CMAKE_C_FLAGS,
   # CMAKE_CXX_FLAGS, etc.
   foreach( comp ${compiler_flag_var_names} )

      # sanity check
      if( NOT ${comp} STREQUAL "C" AND
            NOT ${comp} STREQUAL "CXX" AND
            NOT ${comp} STREQUAL "Fortran" AND
            NOT ${comp} STREQUAL "EXE_LINKER")
         message(FATAL_ERROR "When calling
toggle_compiler_flag(switch, compiler_flag, compiler_flag_var_names),
compiler_flag_var_names must be set to one or more of these valid
names: C;CXX;EXE_LINKER.")
      endif()
      
      string( REPLACE "+" "x" safe_CMAKE_${comp}_FLAGS
         "${CMAKE_${comp}_FLAGS}" )

      if( "${build_modes}x" STREQUAL "x" ) # set flags for all build modes

         if( ${switch} )
            if( NOT "${safe_CMAKE_${comp}_FLAGS}" MATCHES "${safe_compiler_flag}" )
               set( CMAKE_${comp}_FLAGS "${CMAKE_${comp}_FLAGS} ${compiler_flag}" 
                  CACHE STRING "compiler flags" FORCE )
            endif()
         else()
            if( "${safe_CMAKE_${comp}_FLAGS}" MATCHES "${safe_compiler_flag}" )
               string( REPLACE "${compiler_flag}" "" 
                  CMAKE_${comp}_FLAGS ${CMAKE_${comp}_FLAGS} )
               set( CMAKE_${comp}_FLAGS "${CMAKE_${comp}_FLAGS}" 
                  CACHE STRING "compiler flags" FORCE )
            endif()
         endif()

      else() # build_modes listed
         
         foreach( bm ${build_modes} )

            string( REPLACE "+" "x" safe_CMAKE_${comp}_FLAGS_${bm}
               ${CMAKE_${comp}_FLAGS_${bm}} )
            
            if( ${switch} )
               if( NOT "${safe_CMAKE_${comp}_FLAGS_${bm}}" MATCHES 
                     "${safe_compiler_flag}" )
                  set( CMAKE_${comp}_FLAGS_${bm} 
                     "${CMAKE_${comp}_FLAGS_${bm}} ${compiler_flag}" 
                     CACHE STRING "compiler flags" FORCE )
               endif()
            else()
               if( "${safe_CMAKE_${comp}_FLAGS_${bm}}" MATCHES 
                     "${safe_compiler_flag}" )
                  string( REPLACE "${compiler_flag}" "" 
                     CMAKE_${comp}_FLAGS_${bm} ${CMAKE_${comp}_FLAGS_${bm}} )
                  set( CMAKE_${comp}_FLAGS_${bm} 
                     "${CMAKE_${comp}_FLAGS_${bm}}" 
                     CACHE STRING "compiler flags" FORCE )
               endif()
            endif()
            
         endforeach()

      endif()

   endforeach()
endmacro()

#------------------------------------------------------------------------------#
# End config/compiler_env.cmake
#------------------------------------------------------------------------------#
