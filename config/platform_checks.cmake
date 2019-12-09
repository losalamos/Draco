##---------------------------------------------------------------------------##
# file   : platform_checks.cmake
# brief  : Platform Checks for Draco Build System
# note   : Copyright (C) 2016-2019 Triad National Security, LLC.
#          All rights reserved
##---------------------------------------------------------------------------##

include_guard(GLOBAL)

# ----------------------------------------------------------------------------
# Identify machine and save name in ds++/config.h
# ----------------------------------------------------------------------------

function(dbs_set_sitename)

  if( DEFINED SITENAME )
    if( DEFINED SITENAME_FAMILY)
      return()
    else()
      unset(SITENAME)
    endif()
  endif()

  site_name( SITENAME )
  string( REGEX REPLACE "([A-z0-9]+).*" "\\1" SITENAME ${SITENAME} )
  set( SITENAME_FAMILY "unknown" )
  if( ${SITENAME} MATCHES "ba")
    set( SITENAME "Badger" )
    set( SITENAME_FAMILY "CTS-1" )
  elseif( ${SITENAME} MATCHES "ccscs[0-9]+" )
    set( SITENAME_FAMILY "CCS-NET" )
  elseif( ${SITENAME} MATCHES "fi")
    set( SITENAME "Fire" )
    set( SITENAME_FAMILY "CTS-1" )
  elseif( ${SITENAME} MATCHES "ic")
    set( SITENAME "Ice" )
    set( SITENAME_FAMILY "CTS-1" )
  elseif( ${SITENAME} MATCHES "nid")
    if( "$ENV{SLURM_CLUSTER_NAME}" MATCHES "trinity" )
      set( SITENAME "Trinity" )
      set( SITENAME_FAMILY "ATS-1" )
    else()
      set( SITENAME "Trinitite" )
      set( SITENAME_FAMILY "ATS-1" )
    endif()
  elseif( ${SITENAME} MATCHES "sn")
    set( SITENAME "Snow" )
    set( SITENAME_FAMILY "CTS-1" )
  elseif( ${SITENAME} MATCHES "tr")
    set( SITENAME "Trinity" )
      set( SITENAME_FAMILY "ATS-1" )
  elseif( ${SITENAME} MATCHES "tt")
    set( SITENAME "Trinitite" )
      set( SITENAME_FAMILY "ATS-1" )
  endif()
  set( SITENAME ${SITENAME} CACHE "STRING" "Name of the current machine" FORCE)
  set( SITENAME_FAMILY ${SITENAME_FAMILY} CACHE "STRING"
    "Name of the current machine family (ATS-1, CTS-1, etc.)" FORCE)

endfunction()

if( NOT DEFINED CRAY_PE )
  message("
Platform Checks...
")
  dbs_set_sitename()
endif()

#------------------------------------------------------------------------------#
# Sanity checks for Cray Programming Environments
#
# If this is a Cray PE,
# - Set CRAY_PE = TRUE
# - Ensure CMAKE_EXE_LINKER_FLAGS contains "-dynamic"
# - Ensure that the compilers given to cmake are actually Cray compiler
#   wrappers.
#------------------------------------------------------------------------------#
macro( query_craype )

  if( NOT DEFINED CRAY_PE )

    # Is this a Cray machine?
    message( STATUS "Looking to see if we are building in a Cray Environment...")
    if( DEFINED ENV{CRAYPE_VERSION} )
      set( CRAY_PE ON CACHE BOOL
        "Are we building in a Cray Programming Environment?")

      # We expect developers to use the Cray compiler wrappers (especially in
      # setupMPI.cmake). See also
      # https://cmake.org/cmake/help/latest/module/FindMPI.html
      if( NOT "$ENV{CXX}" MATCHES "CC$" OR
          NOT "$ENV{CC}" MATCHES "cc$" OR
          NOT "$ENV{FC}" MATCHES "ftn$" OR
          NOT "$ENV{CRAYPE_LINK_TYPE}" MATCHES "dynamic$" )
        message( FATAL_ERROR
"The build system requires that the Cray compiler wrappers (CC, cc, ftn) be "
" used when configuring this product on a Cray system (CRAY_PE=${CRAY_PE}). The"
" development environment must also support dynamic linking.  The build system "
" thinks you are trying to use:\n"
"  CMAKE_CXX_COMPILER     = ${CMAKE_CXX_COMPILER}\n"
"  CMAKE_C_COMPILER       = ${CMAKE_C_COMPILER}\n"
"  CMAKE_Fortran_COMPILER = ${CMAKE_Fortran_COMPILER}\n"
"  CRAYPE_LINK_TYPE       = $ENV{CRAYPE_LINK_TYPE}\n"
"If you are working on a system that runs the Cray Programming Environment, try"
" setting the following variables and rerunning cmake from a clean build"
" directory:\n"
"   export CXX=`which CC`\n"
"   export CC=`which cc`\n"
"   export FC=`which ftn`\n"
"   export CRAYPE_LINK_TYPE=dynamic\n"
"Otherwise please email this error message and other related information to"
" draco@lanl.gov.\n" )
      endif()
      message( STATUS
        "Looking to see if we are building in a Cray Environment..."
        "found version $ENV{CRAYPE_VERSION}.")
    else()
      set( CRAY_PE OFF CACHE BOOL
        "Are we building in a Cray Programming Environment?")
      message( STATUS
        "Looking to see if we are building in a Cray Environment...no.")
    endif()

  endif()
endmacro()

##---------------------------------------------------------------------------##
## Determine System Type and System Names
##
## Used by ds++ and c4.
##---------------------------------------------------------------------------##
macro( set_draco_uname )
    # Store platform information in config.h
    if( UNIX )
      set( draco_isLinux 1 )
      set( DRACO_UNAME "Linux" )
    elseif( WIN32 )
      set( draco_isWin 1 )
      set( DRACO_UNAME "Windows" )
    elseif( OSF1 )
      set( draco_isOSF1 1 )
      set( DRACO_UNAME "OSF1" )
    elseif( APPLE )
      set( draco_isDarwin 1)
      set( DRACO_UNAME "Darwin" )
    else()
      set( draco_isAIX 1 )
      set( DRACO_UNAME "AIX" )
    endif()
endmacro()

##---------------------------------------------------------------------------##
## Determine if gethostname() is available.
## Determine the value of HOST_NAME_MAX.
##
## Used by ds++/SystemCall.cc and ds++/path.cc
##---------------------------------------------------------------------------##
macro( query_have_gethostname )
    # Platform checks for gethostname()
    include( CheckIncludeFiles )
    check_include_files( unistd.h    HAVE_UNISTD_H )
    check_include_files( limits.h    HAVE_LIMITS_H )
    check_include_files( winsock2.h  HAVE_WINSOCK2_H )
    check_include_files( direct.h    HAVE_DIRECT_H )
    check_include_files( sys/param.h HAVE_SYS_PARAM_H )
    # Used to demangle symbols for stack trace
    # check_include_files( cxxabi.h    HAVE_CXXABI_H )

    # -------------- Checks for hostname and len(hostname) ---------------- #
    # gethostname()
    include( CheckFunctionExists )
    check_function_exists( gethostname HAVE_GETHOSTNAME )

    # HOST_NAME_MAX
    include( CheckSymbolExists )
    unset( hlist )
    if( HAVE_UNISTD_H )
       list( APPEND hlist unistd.h )
    endif()
    if( HAVE_WINSOCK2_H )
       list( APPEND hlist winsock2.h )
    endif()
    if( HAVE_LIMITS_H )
       list( APPEND hlist limits.h )
    endif()
    check_symbol_exists( HOST_NAME_MAX "${hlist}" HAVE_HOST_NAME_MAX )
    if( NOT HAVE_HOST_NAME_MAX )
       unset( HAVE_GETHOSTNAME )
    endif()

    check_symbol_exists( _POSIX_HOST_NAME_MAX "posix1_lim.h"
      HAVE_POSIX_HOST_NAME_MAX )

    # HOST_NAME_MAX
    check_symbol_exists( MAXHOSTNAMELEN "sys/param.h" HAVE_MAXHOSTNAMELEN )
    if( NOT HAVE_MAXHOSTNAMELEN )
       unset( HAVE_MAXHOSTNAMELEN )
    endif()

endmacro()

##---------------------------------------------------------------------------##
## Determine if gethostname() is available.
## Determine the value of HOST_NAME_MAX.
##
## Used by ds++/SystemCall.cc and ds++/path.cc
##---------------------------------------------------------------------------##
macro( query_have_maxpathlen )
    # MAXPATHLEN
    unset( hlist )
    if( HAVE_UNISTD_H )
       list( APPEND hlist unistd.h )
    endif()
    if( HAVE_LIMITS_H )
       list( APPEND hlist limits.h )
    endif()
    if( HAVE_SYS_PARAM_H )
       list( APPEND hlist sys/param.h )
    endif()
    check_symbol_exists( MAXPATHLEN "${hlist}" HAVE_MAXPATHLEN )
    if( NOT HAVE_MAXPATHLEN )
        unset( HAVE_MAXPATHLEN )
    endif()
endmacro()

##---------------------------------------------------------------------------##
## Determine if some system headers exist
##---------------------------------------------------------------------------##
macro( query_have_sys_headers )

   include( CheckIncludeFiles )
   check_include_files( sys/types.h HAVE_SYS_TYPES_H )
   check_include_files( unistd.h    HAVE_UNISTD_H    )

endmacro()

##---------------------------------------------------------------------------##
## Check 8-byte int type
##
## For some systems, provide special compile flags to support 8-byte integers
##---------------------------------------------------------------------------##

macro(check_eight_byte_int_type)
   if( "${SIZEOF_INT}notset" STREQUAL "notset" )
      determine_word_types()
   endif()

   if( "${SIZEOF_INT}" STREQUAL "8" )
      message( "Checking for 8-byte integer type... int - no mods needed." )
   elseif( "${SIZEOF_LONG}" STREQUAL "8" )
      message( "Checking for 8-byte integer type... long - no mods needed." )
   else()
      message( FATAL_ERROR "need to patch up this part of the build system." )
   endif()
endmacro()

##---------------------------------------------------------------------------##
## Detect support for the C99 restrict keyword
## Borrowed from http://cmake.3232098.n2.nabble.com/AC-C-RESTRICT-td7582761.html
##
## A restrict-qualified pointer (or reference) is basically a promise to the
## compiler that for the scope of the pointer, the target of the pointer will
## only be accessed through that pointer (and pointers copied from it).
##
## http://www.research.scea.com/research/pdfs/GDC2003_Memory_Optimization_18Mar03.pdf
##---------------------------------------------------------------------------##
macro( query_have_restrict_keyword )

  if( NOT PLATFORM_CHECK_RESTRICT_KYEWORD_DONE )
    set( PLATFORM_CHECK_RESTRICT_KYEWORD_DONE TRUE CACHE BOOL
      "Is restrict keyword check done?")
    mark_as_advanced( PLATFORM_CHECK_RESTRICT_KYEWORD_DONE )
    message(STATUS "Looking for the C99 restrict keyword")
    include( CheckCSourceCompiles )
    foreach( ac_kw __restrict __restrict__ _Restrict restrict )
      check_c_source_compiles("
         typedef int * int_ptr;
         int foo ( int_ptr ${ac_kw} ip ) { return ip[0]; }
         int main (void) {
            int s[1];
            int * ${ac_kw} t = s;
            t[0] = 0;
            return foo(t); }
         "
        HAVE_RESTRICT)

      if( HAVE_RESTRICT )
        set( RESTRICT_KEYWORD ${ac_kw} )
        message(STATUS
          "Looking for the C99 restrict keyword - found ${RESTRICT_KEYWORD}")
        break()
      endif()
    endforeach()
    if( NOT HAVE_RESTRICT )
      message(STATUS "Looking for the C99 restrict keyword - not found")
    endif()
  endif()

endmacro()

#------------------------------------------------------------------------------#
# Query if hardware has FMA, AVX2
#
# This code is adopted from
# https://software.intel.com/en-us/node/405250?language=es&wapkw=avx2+cpuid
#------------------------------------------------------------------------------#
macro( query_fma_on_hardware )

  if( NOT PLATFORM_CHECK_FMA_DONE )

    set( PLATFORM_CHECK_FMA_DONE TRUE CACHE BOOL
      "Is the check for hardware FMA done?")
    mark_as_advanced( PLATFORM_CHECK_FMA_DONE )
    message( STATUS "Looking for hardware FMA support...")

    if( "${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "ppc64le" OR
        "${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "aarch64" )
      # recent arm and power8/9 chips have FMA and the check below fails for
      # these architectures, so we hard code the result here.
      set(HAVE_HARDWARE_FMA TRUE)

    else()
      unset(HAVE_HARDWARE_FMA)
      try_run(
        HAVE_HARDWARE_FMA
        HAVE_HARDWARE_FMA_COMPILE
        ${CMAKE_CURRENT_BINARY_DIR}/config
        ${CMAKE_CURRENT_SOURCE_DIR}/config/query_fma.cc
        ARGS "-f"
        )
      if( NOT HAVE_HARDWARE_FMA_COMPILE )
        message( FATAL_ERROR "Unable to compile config/query_fma.cc.")
      endif()
    endif()

    if( HAVE_HARDWARE_FMA )
      message( STATUS "Looking for hardware FMA support...found fma.")
    else()
      message( STATUS "Looking for hardware FMA support...fma not found.")
    endif()

    # Other things to look at (might be able to avoid the try-compile):

    # if (WIN32)
    #   # Not sure what to do here. Consider:
    #   # - looking at $ENV{PROCESSOR_IDENTIFIER}. This will be something like:
    #   #   "Intel64 Family 6 Model 45 Stepping 7, GenuineIntel" This string
    #   #   would need to be decoded to know if the processor supports FMA.
    #   # - running 'wmic cpu get * /fomrat:list'. This lists a lot of
    #   #   information about the cpu, but it does not itemize features like
    #   #   fma. Optionally, 'wmic cpu get name'
    #   # - run a 3rd party application like cpuz64.
    # elseif (APPLE)
    #   execute_process( COMMAND /usr/sbin/sysctl -n hw.optional.fma
    #     OUTPUT_VARIABLE found_fma
    #     OUTPUT_QUIET )
    #   if( ${found_fma} GREATER 0 )
    #     set(HAS_HARDWARE_FMA ON)
    #   endif()
    # else()
    #   if( EXISTS /proc/cpuinfo )
    #     execute_process( COMMAND /bin/cat /proc/cpuinfo
    #       OUTPUT_VARIABLE cpuinfo-output )
    #   string( FIND "${cpuinfo-output}" fma found_fma )
    #   if( ${found_fma} GREATER 0 )
    #     set(HAS_HARDWARE_FMA ON)
    #   endif()
    # endif()

    message( STATUS "Looking for hardware AVX2 support...")

    if( "${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "ppc64le" )
      # see comments above for FMA
      set(HAVE_HARDWARE_AVX2 TRUE)

    elseif( "${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "aarch64")
      # see comments above for FMA
      set(HAVE_HARDWARE_AVX2 FALSE)

    else()
      unset(HAVE_HARDWARE_AVX2)
      try_run(
        HAVE_HARDWARE_AVX2
        HAVE_HARDWARE_AVX2_COMPILE
        ${CMAKE_CURRENT_BINARY_DIR}/config
        ${CMAKE_CURRENT_SOURCE_DIR}/config/query_fma.cc
        ARGS "-f"
        )
      if( NOT HAVE_HARDWARE_AVX2_COMPILE )
        message( FATAL_ERROR "Unable to compile config/query_fma.cc.")
      endif()
    endif()

    if( HAVE_HARDWARE_AVX2 )
      message( STATUS "Looking for hardware AVX2 support...found AVX2.")
    else()
      message( STATUS "Looking for hardware AVX2 support...AVX2 not found.")
    endif()

  endif()

endmacro()

##---------------------------------------------------------------------------##
## Sample platform checks
##---------------------------------------------------------------------------##

# # Check for nonblocking collectives
# check_function_exists(MPI_Iallgather HAVE_MPI3_NONBLOCKING_COLLECTIVES)
# check_function_exists(MPIX_Iallgather HAVE_MPIX_NONBLOCKING_COLLECTIVES)

# # Check for MPI_IN_PLACE (essentially MPI2 support)
# include(CheckCXXSourceCompiles)
# set(MPI_IN_PLACE_CODE
#     "#include \"mpi.h\"
#      int main( int argc, char* argv[] )
#      {
#          MPI_Init( &argc, &argv );
#          float a;
#          MPI_Allreduce
#          ( MPI_IN_PLACE, &a, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD );
#          MPI_Finalize();
#          return 0;
#      }
#     ")
# set(CMAKE_REQUIRED_FLAGS ${CXX_FLAGS})
# set(CMAKE_REQUIRED_INCLUDES ${MPI_CXX_INCLUDE_PATH})
# set(CMAKE_REQUIRED_LIBRARIES ${MPI_CXX_LIBRARIES})
# check_cxx_source_compiles("${MPI_IN_PLACE_CODE}" HAVE_MPI_IN_PLACE)

# # Look for restrict support
# set(RESTRICT_CODE
#     "int main(void)
#      {
#          int* RESTRICT a;
#          return 0;
#      }")
# set(CMAKE_REQUIRED_DEFINITIONS "-DRESTRICT=__restrict__")
# check_cxx_source_compiles("${RESTRICT_CODE}" HAVE___restrict__)
# if(HAVE___restrict__)
# ...
# endif()



##---------------------------------------------------------------------------##
