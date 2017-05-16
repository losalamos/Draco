##---------------------------------------------------------------------------##
# file   : platform_checks.cmake
# brief  : Platform Checks for Draco Build System
# note   : Copyright (C) 2016-2017 Los Alamos National Security, LLC.
#          All rights reserved
##---------------------------------------------------------------------------##

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

    check_symbol_exists( _POSIX_HOST_NAME_MAX "posix1_lim.h" HAVE_POSIX_HOST_NAME_MAX )

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
         message(STATUS "Looking for the C99 restrict keyword - found ${RESTRICT_KEYWORD}")
         break()
      endif()
   endforeach()
   if( NOT HAVE_RESTRICT )
      message(STATUS "Looking for the C99 restrict keyword - not found")
   endif()

endmacro()

##---------------------------------------------------------------------------##
## Query OpenMP availability
##
## This feature is usually compiler specific and a compile flag must be
## added.  For this to work the <platform>-<compiler>.cmake files (eg.
## unix-g++.cmake) call this macro.
##---------------------------------------------------------------------------##
macro( query_openmp_availability )
  message( STATUS "Looking for OpenMP...")
  find_package(OpenMP QUIET)
  if( OPENMP_FOUND )
    message( STATUS "Looking for OpenMP... ${OpenMP_C_FLAGS}")
    set( OPENMP_FOUND ${OPENMP_FOUND} CACHE BOOL "Is OpenMP availalable?" FORCE )
  else()
    message(STATUS "Looking for OpenMP... not found")
  endif()
endmacro()

#------------------------------------------------------------------------------#
# Query if hardware has FMA
#
# This code is adopted from
# https://software.intel.com/en-us/node/405250?language=es&wapkw=avx2+cpuid
#------------------------------------------------------------------------------#
macro( query_fma_on_hardware )

  message( STATUS "Looking for hardware FMA support...")
  unset(HAVE_HARDWARE_FMA)
  try_run(
    HAVE_HARDWARE_FMA
    HAVE_HARDWARE_FMA_COMPILE
    ${CMAKE_CURRENT_BINARY_DIR}/config
    ${CMAKE_CURRENT_SOURCE_DIR}/config/query_fma.cc
    )
  if( NOT HAVE_HARDWARE_FMA_COMPILE )
    message( FATAL_ERROR "Unable to compile config/query_fma.cc.")
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
  #   #   "Intel64 Family 6 Model 45 Stepping 7, GenuineIntel" This string would
  #   #   need to be decoded to know if the processor supports FMA.
  #   # - running 'wmic cpu get * /fomrat:list'. This lists a lot of information
  #   #   about the cpu, but it does not itemize features like fma. Optionally,
  #   #   'wmic cpu get name'
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
