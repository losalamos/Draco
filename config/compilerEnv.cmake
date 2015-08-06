#-----------------------------*-cmake-*----------------------------------------#
# file   config/compilerEnv.cmake
# brief  Default CMake build parameters
# note   Copyright (C) 2010-2015 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

include( FeatureSummary )

# ----------------------------------------
# PAPI
# ----------------------------------------
if( EXISTS $ENV{PAPI_HOME} )

  set( HAVE_PAPI 1 CACHE BOOL "Is PAPI available on this machine?" )
  set( PAPI_INCLUDE $ENV{PAPI_INCLUDE} CACHE PATH "PAPI headers at this location" )
  set( PAPI_LIBRARY $ENV{PAPI_LIBDIR}/libpapi.so CACHE FILEPATH "PAPI library." )
endif()

# PAPI 4.2 on CT uses a different setup.
if( $ENV{PAPI_VERSION} MATCHES "[45].[0-9].[0-9]")
  set( HAVE_PAPI 1 CACHE BOOL "Is PAPI available on this machine?" )
  string( REGEX REPLACE ".*[ ][-]I(.*)$" "\\1" PAPI_INCLUDE $ENV{PAPI_INCLUDE_OPTS} )
  string( REGEX REPLACE ".*[ ][-]L(.*)[ ].*" "\\1" PAPI_LIBDIR $ENV{PAPI_POST_LINK_OPTS} )
endif()

if( HAVE_PAPI )
  set( PAPI_INCLUDE ${PAPI_INCLUDE} CACHE PATH "PAPI headers at this location" )
  set( PAPI_LIBRARY ${PAPI_LIBDIR}/libpapi.so CACHE FILEPATH "PAPI library." )
  if( NOT EXISTS ${PAPI_LIBRARY} )
    message( FATAL_ERROR "PAPI requested, but library not found.  Set PAPI_LIBDIR to correct path." )
  endif()
  mark_as_advanced( PAPI_INCLUDE PAPI_LIBRARY )
  add_feature_info( HAVE_PAPI HAVE_PAPI "Provide PAPI hardware counters if available." )
endif()

# ----------------------------------------
# MIC processors
# ----------------------------------------
if( "${HAVE_MIC}x" STREQUAL "x" )

  # default to OFF
  set( HAVE_MIC OFF)

  # This was the old mechanism.  It fails to work because we might be
  # targeting a haswell node that also has MIC processors. Still, we
  # might want to use this in the future to determine if MICs are on
  # the local node.
  message( STATUS "Looking for availability of MIC hardware")
  set( mic_found FALSE )
  if( EXISTS /usr/sbin/micctrl )
    exec_program(
      /usr/sbin/micctrl
      ARGS -s
      OUTPUT_VARIABLE mic_status
      OUTPUT_STRIP_TRAILING_WHITESPACE
      OUTPUT_QUIET
      )
    if( mic_status MATCHES online )
      # set( HAVE_MIC ON CACHE BOOL "Does the local machine have MIC chips?" )
      # If we areusing MIC, then disable CUDA.
      # set( USE_CUDA OFF CACHE BOOL "Compile against Cuda libraries?")
      set( mic_found TRUE )
    endif()
  endif()
  if( mic_found )
    message( STATUS "Looking for availability of MIC hardware - found")
  else()
    message( STATUS "Looking for availability of MIC hardware - not found")
  endif()

  if( ${mic_found} )
    # Should we cross compile for the MIC?
    # Look at the environment variable SLURM_JOB_PARTITION to determine
    # if we should cross compile for the MIC processor.
    message(STATUS "Enable cross compiling for MIC (HAVE_MIC) ...")
    # See https://darwin.lanl.gov/darwin_hw/report.html for a list
    # of string designators used as partition names:
    if( NOT "$ENV{SLURM_JOB_PARTITION}x" STREQUAL "x" AND "$ENV{SLURM_JOB_PARTITION}" STREQUAL "knc-mic")
      set( HAVE_MIC ON )
      set( USE_CUDA OFF CACHE BOOL "Compile against Cuda libraries?")
    endif()

    # Store the result in the cache.
    set( HAVE_MIC ${HAVE_MIC} CACHE BOOL "Should we cross compile for the MIC processor?" )
    message(STATUS "Enable cross compiling for MIC (HAVE_MIC) ... ${HAVE_MIC}")
    unset( mic_found)
  endif()

endif()

# ----------------------------------------
# STATIC or SHARED libraries?
# ----------------------------------------

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

endmacro()

#------------------------------------------------------------------------------#
# Setup C++ Compiler
#------------------------------------------------------------------------------#
macro(dbsSetupCxx)

  dbsSetupCompilers()

  # Deal with compiler wrappers
  if( ${CMAKE_CXX_COMPILER} MATCHES "tau_cxx.sh" )
    # When using the TAU profiling tool, the actual compiler vendor
    # is hidden under the tau_cxx.sh script.  Use the following
    # command to determine the actual compiler flavor before setting
    # compiler flags (end of this macro).
    execute_process(
      COMMAND ${CMAKE_CXX_COMPILER} -tau:showcompiler
      OUTPUT_VARIABLE my_cxx_compiler )
  elseif( ${CMAKE_CXX_COMPILER} MATCHES "xt-asyncpe" OR
          ${CMAKE_CXX_COMPILER} MATCHES "craype" )
    # Ceilo (catamount) uses a wrapper script
    # /opt/cray/xt-asyncpe/5.06/bin/CC that masks the actual compiler.
    # Use the following command to determine the actual compiler
    # flavor before setting compiler flags (end of this macro).
    # Trinitite uses a similar wrapper script (/opt/cray/craype/2.4.0/bin/CC).
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

  # C99 support:
  set( CMAKE_C_STANDARD 99 )

  # C++11 support:
  set( CMAKE_CXX_STANDARD 11 )

  # Do not enable extensions (e.g.: --std=gnu++11)
  set( CMAKE_CXX_EXTENSIONS OFF )
  set( CMAKE_C_EXTENSIONS   OFF )

  get_filename_component( my_cxx_compiler ${my_cxx_compiler} NAME )
  if( ${my_cxx_compiler} MATCHES "mpicxx" )
    # MPI wrapper
    execute_process( COMMAND ${my_cxx_compiler} --version
      OUTPUT_VARIABLE mpicxx_version_output
      OUTPUT_STRIP_TRAILING_WHITESPACE )
    # make output safe for regex matching
    string( REPLACE "+" "x" mpicxx_version_output ${mpicxx_version_output} )
    if( ${mpicxx_version_output} MATCHES icpc )
      set( my_cxx_compiler icpc )
    endif()
  endif()
  if( ${my_cxx_compiler} MATCHES "clang" OR
      ${my_cxx_compiler} MATCHES "llvm")
    if( APPLE )
      include( apple-clang )
    else()
      include( unix-clang )
    endif()
  elseif( ${my_cxx_compiler} MATCHES "cl" )
    include( windows-cl )
  elseif( ${my_cxx_compiler} MATCHES "icpc" )
    include( unix-intel )
  elseif( ${my_cxx_compiler} MATCHES "pgCC" )
    include( unix-pgi )
  elseif( ${my_cxx_compiler} MATCHES "xl" )
    include( unix-xl )
  elseif( ${my_cxx_compiler} MATCHES "[cg][+x]+" )
    include( unix-g++ )
  else()
    message( FATAL_ERROR "Build system does not support CXX=${my_cxx_compiler}" )
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
  add_definitions(-D_POSIX_C_SOURCE=200112 -D_XOPEN_SOURCE=600)
  set( CMAKE_REQUIRED_DEFINITIONS "${CMAKE_REQUIRED_DEFINITIONS} -D_POSIX_C_SOURCE=200112" )
  set( CMAKE_REQUIRED_DEFINITIONS "${CMAKE_REQUIRED_DEFINITIONS} -D_XOPEN_SOURCE=600")
  if ( APPLE )
    # Defining the above requires adding POSIX extensions,
    # otherwise, include ordering still goes wrong on Darwin,
    # (i.e., putting fstream before iostream causes problems)
    # see https://code.google.com/p/wmii/issues/detail?id=89
    add_definitions(-D_DARWIN_C_SOURCE)
    set( CMAKE_REQUIRED_DEFINITIONS "${CMAKE_REQUIRED_DEFINITIONS} -D_DARWIN_C_SOURCE ")
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

  dbsSetupCompilers()

  # Toggle if we should try to build Fortran parts of the project.
  # This will be set to true if $ENV{FC} points to a working compiler
  # (e.g.: GNU or Intel compilers with Unix Makefiles) or if the
  # current project doesn't support Fortran but
  # CMakeAddFortranSubdirectory can be used.
  option( HAVE_Fortran "Should we build Fortran parts of the project?" OFF )

  # Is Fortran enabled (it is considered 'optional' for draco)?
  get_property(_LANGUAGES_ GLOBAL PROPERTY ENABLED_LANGUAGES)
  if( _LANGUAGES_ MATCHES Fortran )

    set( HAVE_Fortran ON )

    # Deal with comiler wrappers
    if( ${CMAKE_Fortran_COMPILER} MATCHES "xt-asyncpe" OR
        ${CMAKE_Fortran_COMPILER} MATCHES "craype" )
      # Ceilo (catamount) uses a wrapper script
      # /opt/cray/xt-asyncpe/5.06/bin/CC that masks the actual
      # compiler.  Use the following command to determine the actual
      # compiler flavor before setting compiler flags (end of this
      # macro).
      # Trinitite uses a similar wrapper script (/opt/cray/craype/2.4.0/bin/ftn).
      execute_process(
        COMMAND ${CMAKE_Fortran_COMPILER} --version
        OUTPUT_VARIABLE my_fc_compiler
        ERROR_QUIET )
      string( REGEX REPLACE "^(.*).Copyright.*" "\\1"
        my_fc_compiler ${my_fc_compiler})
    else()
      set( my_fc_compiler ${CMAKE_Fortran_COMPILER} )
    endif()

    # MPI wrapper
    if( ${my_fc_compiler} MATCHES "mpif90" )
      execute_process( COMMAND ${my_fc_compiler} --version
        OUTPUT_VARIABLE mpifc_version_output
        OUTPUT_STRIP_TRAILING_WHITESPACE )
      if( ${mpifc_version_output} MATCHES ifort )
        set( my_fc_compiler ifort )
      elseif( ${mpifc_version_output} MATCHES GNU )
        set( my_fc_compiler gfortran )
      endif()
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
      message( FATAL_ERROR "Build system does not support FC=${my_fc_compiler}" )
    endif()

    if( _LANGUAGES_ MATCHES "^C$" OR _LANGUAGES_ MATCHES CXX )
      include(FortranCInterface)
    endif()

  else()
    # If CMake doesn't know about a Fortran compiler, $ENV{FC}, then
    # also look for a compiler to use with
    # CMakeAddFortranSubdirectory.
    message( STATUS "Looking for CMakeAddFortranSubdirectory Fortran compiler...")

    # Try to find a Fortran compiler (use MinGW gfortran for MSVC).
    find_program( CAFS_Fortran_COMPILER
      NAMES ${CAFS_Fortran_COMPILER} $ENV{CAFS_Fortran_COMPILER} gfortran
      PATHS
      c:/MinGW/bin
      "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\Windows\\CurrentVersion\\Uninstall\\MinGW;InstallLocation]/bin" )

    if( EXISTS ${CAFS_Fortran_COMPILER} )
      set( HAVE_Fortran ON )
      message( STATUS "Looking for CMakeAddFortranSubdirectory Fortran compiler... found ${CAFS_Fortran_COMPILER}")
    else()
      message( STATUS "Looking for CMakeAddFortranSubdirectory Fortran compiler... not found")
    endif()

  endif()

endmacro()

##---------------------------------------------------------------------------##
## Toggle a compiler flag based on a bool
##
## Examples:
##   toggle_compiler_flag( GCC_ENABLE_ALL_WARNINGS "-Weffc++" "CXX" "DEBUG" )
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
          set( CMAKE_${comp}_FLAGS "${CMAKE_${comp}_FLAGS} ${compiler_flag} "
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
              "${CMAKE_${comp}_FLAGS_${bm}} ${compiler_flag} "
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
