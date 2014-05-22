#-----------------------------*-cmake-*----------------------------------------#
# file   windows-cl.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 June 5
# brief  Establish flags for Windows - MSVC
# note   Copyright (C) 2010-2014 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$ 
#------------------------------------------------------------------------------#

# Prepare
string( TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE )

# OpenMP is not available in free MSVC products.
if( USE_OPENMP )
    # Platform checks for gethostname()
    include( CheckIncludeFiles )
    check_include_files( omp.h HAVE_OMP_H )
    if( NOT HAVE_OMP_H )
       set( USE_OPENMP OFF CACHE BOOL "Turn on OpenMP features?" FORCE )
    endif()
endif()

# This is required to provide compatibility between MSVC and MinGW generated libraries.
if( DRACO_SHARED_LIBS )
  set( CMAKE_GNUtoMS ON CACHE BOOL "Compatibility flag for MinGW/MSVC." FORCE)
endif()

# Extra setup (ds++/config.h) for MSVC
# 1. Allow M_PI to be found via <cmath>
set( _USE_MATH_DEFINES 1 )

  # MSVC 9 2008 SP1 Flags
  # http://msdn.microsoft.com/en-us/library/19z1t1wy.aspx
  # /arch:SSE2 - enable use of SSE2 instructions (not compatible with /fp:strict)
  # /EHsc - The exception-handling model that catches C++ exceptions only and
  #       tells the compiler to assume that extern C functions never throw 
  #       a C++ exception.
  # /EHa - The exception-handling model that catches asynchronous (structured)
  #       and synchronous (C++) exceptions.
  # /fp:strict - The strictest floating-point model.
  # /GR - Enables run-time type information (RTTI).
  # /Gy - Enables function-level linking
  # /GZ - Enable stack checks (deprecated, use /RTC1 instead, http://msdn.microsoft.com/en-us/library/hddybs7t.aspx)
  # /Od - Disable optimization.
  # /O2 - Create fast code.
  # /01 - Create small code.
  # /Ob - Controls inline expansion
  # /openmp - Enables #pragma omp in source code.
  # /MD - Creates a multithreaded DLL using MSVCRT.lib
  # /MP<N> - Enable parallel builds
  # /MT - Creates a multithreaded executable file using LIBCMT.lib
  # /RTC1 - Enable run-time error checking.
  # /W4 - Issue all warnings.
  # /Za - Disables language extensions. Emits an error for language constructs that are not compatible with either ANSI C or ANSI C++. (DOM parser fails to compile with this flag).
  # /Zi - Generates complete debugging information.

  set( MD_or_MT_debug "${MD_or_MT}d" )
  if( "${DEBUG_RUNTIME_EXT}" STREQUAL "d" )
    set( MD_or_MT_debug "${MD_or_MT}${DEBUG_RUNTIME_EXT} /RTC1" ) # RTC requires /MDd
  endif()

  set( numproc $ENV{NUMBER_OF_PROCESSORS} )
  if( "${numproc}notfound" STREQUAL "notfound" )
    set( numproc 1 )
  endif( "${numproc}notfound" STREQUAL "notfound" )     
  
if( ${CMAKE_GENERATOR} STREQUAL "Visual Studio 8 2005" OR 
   ${CMAKE_GENERATOR}  STREQUAL "Visual Studio 9 2008" OR 
   ${CMAKE_GENERATOR}  STREQUAL "Visual Studio 10" OR
   ${CMAKE_GENERATOR}  STREQUAL "Visual Studio 11" OR
   ${CMAKE_GENERATOR}  STREQUAL "Visual Studio 12" OR
   ${CMAKE_GENERATOR}  STREQUAL "Visual Studio 12 2013" OR
   ${CMAKE_GENERATOR}  STREQUAL "Visual Studio 11 2012" OR
   ${CMAKE_GENERATOR}  MATCHES  "NMake Makefiles" )
   
  set( CMAKE_C_FLAGS "/W2 /Gy /DWIN32 /D_WINDOWS /MP${numproc}" ) # /Za
  set( CMAKE_C_FLAGS_DEBUG "/${MD_or_MT_debug} /Od /Zi /Ob0 /DDEBUG /D_DEBUG" )
  set( CMAKE_C_FLAGS_RELEASE "/${MD_or_MT} /O2 /Ob2 /DNDEBUG" ) # /O2 /Ob2
  set( CMAKE_C_FLAGS_MINSIZEREL "/${MD_or_MT} /O1 /Ob1 /DNDEBUG" )
  set( CMAKE_C_FLAGS_RELWITHDEBINFO "/${MD_or_MT} /O2 /Ob2 /Zi /DDEBUG" ) # /O2 /Ob2

  set( CMAKE_CXX_FLAGS "/W2 /Gy /EHa /DWIN32 /D_WINDOWS /MP${numproc}" ) # /Zm1000 /GX /GR /Za
  set( CMAKE_CXX_FLAGS_DEBUG "/${MD_or_MT_debug} /Od /Zi /Ob0 /DDEBUG /D_DEBUG" )
  set( CMAKE_CXX_FLAGS_RELEASE "/${MD_or_MT} /O2 /Ob2 /DNDEBUG" ) # /O2 /Ob2
  set( CMAKE_CXX_FLAGS_MINSIZEREL "/${MD_or_MT} /O1 /Ob1 /DNDEBUG" )
  set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "/${MD_or_MT} /O2 /Ob2 /Zi /DDEBUG" ) # /O2 /Ob2
  
  # If building static libraries, inlcude debugging information in the
  # library.
  if( ${DRACO_LIBRARY_TYPE} MATCHES "STATIC" )
    set( CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS_DEBUG} /Z7"   )
    set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /Z7" )
  endif()  

  # Suppress some MSVC warnings about "unsafe" pointer use.
  if(MSVC_VERSION GREATER 1399)
    set( CMAKE_C_FLAGS 
      "${CMAKE_C_FLAGS} /D_CRT_SECURE_NO_DEPRECATE /D_SCL_SECURE_NO_DEPRECATE /D_SECURE_SCL=0" )
    set( CMAKE_CXX_FLAGS 
      "${CMAKE_CXX_FLAGS} /D_CRT_SECURE_NO_DEPRECATE /D_SCL_SECURE_NO_DEPRECATE /D_SECURE_SCL=0" )
    #set( CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS_DEBUG} /D_HAS_ITERATOR_DEBUGGING=0" )
    #set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /D_HAS_ITERATOR_DEBUGGING=0" )
  endif(MSVC_VERSION GREATER 1399)

  find_library( Lib_win_winsock 
     NAMES wsock32;winsock32;ws2_32 
     HINTS 
        "C:/Program Files (x86)/Microsoft SDKs/Windows/v7.0A/Lib"
        "C:/Windows/SysWOW64"
  )  
  if( ${Lib_win_winsock} MATCHES "NOTFOUND" )
     message( FATAL_ERROR "Could not find library winsock32 or ws2_32!" )
  endif()
  
else() 

  message( FATAL_ERROR "config/windows-cl.cmake must be taught to build for this compiler (CMAKE_GENERATOR = ${CMAKE_GENERATOR}).  Yell at kt for help on this error." )
  
endif()

#------------------------------------------------------------------------------#
# End windows-cl.cmake
#------------------------------------------------------------------------------#
