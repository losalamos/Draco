#-----------------------------*-cmake-*----------------------------------------#
# file   config/windows-ifort.cmake
# author Kelly Thompson
# date   2008 May 30
# brief  Establish flags for Windows - Intel Visual Fortran
# note   © Copyright 2010 Los Alamos National Security, LLC, All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

# Let anyone who is interested in which FORTRAN compiler we're using
# switch on this macro.
set( CMAKE_Fortran_COMPILER_FLAVOR "IFORT" )

option( ENABLE_Fortran_CODECOVERAGE "Instrument for code coverage analysis?" OFF )
if( ENABLE_Fortran_CODECOVERAGE )
  if( NOT PROF_DIR )
    message( STATUS "Setting PROF_DIR = ${CMAKE_INSTALL_PREFIX}/codecov" )
    set( PROF_DIR ${CMAKE_INSTALL_PREFIX}/codecov )
    set( ENV{PROF_DIR} ${CMAKE_INSTALL_PREFIX}/codecov )
  endif( NOT PROF_DIR )
endif( ENABLE_Fortran_CODECOVERAGE )

   # Suggested by IVF/MSVC defaults
   # 
   # /W1 - enable all warnings
   # /free == free format FORTRAN
   # /fpp == preprocess this file.
   # /G7 Optimizes for Intel Core Duo processors
   # /extfpp:<ext> - Consider file extension <ext> as needing preprocessing.   
   # /fpscomp:logical - Use 1 as .true. (instead of -1) for compatibility
   #              with G95.
   # /fp:strict - Enables value-safe optimizations on floating-point data 
   #              and rounds intermediate results to source-defined precision.
   # /Qprec-div - Attempts to use slower but more accurate implementation of 
   #              floating-point divide. 
   # /Qax       - Tells the compiler to generate multiple, processor-specific
   #              auto-dispatch code paths for Intel processors if there is 
   #              a performance benefit
   # /arch:SSE2 - Generates code for Intel Streaming SIMD Extensions 2.
   # /QaxSSSE3  - Tells the compiler to generate multiple, processor-specific
   #              auto-dispatch code paths for Intel processors if there is a
   #              performance benefit.  Use SSSE3 if it is available.
   # /iface:stdcall - Use STDCALL conventions (@FOO_NAME)
   set( CMAKE_Fortran_FLAGS "/warn /free /fpp /extfpp:.F95 /fpscomp:logical /fp:strict" ) # 2009-12-31 (KT) /fp:stric is important for passing the fireball tests.
   # Before 2009-11-04: "/W1 /free /fpp /G7 /extfpp:.F95 /fpscomp:logical /fp:strict"
   if( MSVC_IDE )
     set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /nologo" )
   endif( MSVC_IDE )

   # /MDd == resolve system calls with multi-threaded, dynamic link libraries (debug)
   # /debug:full (==/Zi)
   # /check:bounds - 
   # /check:all
   # /traceback - Tells the compiler to generate extra information in
   #              the object file to provide source file traceback
   #              information when a severe error occurs at run time.
   # /Qtrapuv   - Trap uninitialized variables.
   # /Quse-vcdebug - Issue debug information that is compatible with
   #              MSVC debugger (deprecated after version 10.1)
   # /Qvec-report0 (don't report on vectorization diagnostics/progress.)
   set(CMAKE_Fortran_FLAGS_DEBUG 
     "/Od /debug:full  /traceback /Qtrapuv /DDEBUG /${MD_or_MT}${DEBUG_RUNTIME_EXT}" )  # /Quse-vcdebug 
   set(CMAKE_Fortran_FLAGS_RELEASE "/O2 /Ob2 /fp:except- /Qvec-report:0 /Qopenmp-report:0 /Qpar-report:0 /DNDEBUG /${MD_or_MT}" ) # /${MD_or_MT} /O3
   set(CMAKE_Fortran_FLAGS_MINSIZEREL "/Os /fp:except- /Qvec-report:0 /Qopenmp-report:0 /Qpar-report:0 /DNDEBUG /${MD_or_MT}" )  # /${MD_or_MT}
   set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "/O2 /Ob2 /fp:except- /DDEBUG /${MD_or_MT}" )  # /${MD_or_MT} /Quse-vcdebug  /traceback /debug:extended /debug:full
   # Before 2009-11-04: "/O2 /Ob2 /debug:full /traceback /debug:extended /DDEBUG /${MD_or_MT}"
   
   # Save the compiler version
   execute_process( 
      COMMAND ${CMAKE_Fortran_COMPILER}
      ERROR_VARIABLE tmp ) 
   string( REGEX REPLACE "\n" " " tmp "${tmp}" )
   string( REGEX REPLACE "(.*ID[: ]+)([bcfoprw_]+)([0-9]+[.][0-9]+[.][0-9]+)(.*)" "\\3"
      CMAKE_Fortran_COMPILER_VERSION "${tmp}" )
     
   # save the major version number
   string( REGEX REPLACE "([0-9]+).*" "\\1" CMAKE_Fortran_COMPILER_VERSION_MAJOR 
     "${CMAKE_Fortran_COMPILER_VERSION}" )
 
#
# During discovery of F95 compiler, also discover and make available:
#
# ${CMAKE_Fortran_redist_dll}    
#   - List of Fortran compiler libraries to be installed for release
# ${CMAKE_Fortran_debug_dll}
#   - List of Fortran compiler libraries to be installed for portable developer-only debug version
# ${CMAKE_Fortran_compiler_libs} 
#   - List of Fortran compiler libraries to be used with the target_link_libraries command (C main code that links with Fortran built library.)

# ONLY non-debug versions are redistributable.
set( f90_system_dll 
   libifcoremd.dll
   libifportmd.dll
   libmmd.dll
   )
   
set( f90_system_lib
  #ifconsol.lib
  #libguide40.dll
  #libifcoremd.lib # libifcoremdd.lib
  #libifportmd.lib
  #libirc.lib
  #libmmdd.lib     # libmmd.lib
  #svml_disp.lib   # svml_dispmd.lib?
  )

# Add the correct OpenMP library: older/newer versions, static/dynamic versions, performance/profiling versions:
# (ref: http://www.intel.com/software/products/compilers/docs/flin/main_for/mergedprojects/optaps_for/common/optaps_par_libs.htm)
if( ${CMAKE_Fortran_COMPILER_VERSION_MAJOR} GREATER 10 )
  if( ${MD_or_MT} MATCHES "MD" )
    if( ENABLE_THREAD_PROFILE )
      list( APPEND f90_system_dll libiompprof5md.dll )
      list( APPEND f90_system_lib libiompprof5md.lib )
    else( ENABLE_THREAD_PROFILE )
      list( APPEND f90_system_dll libiomp5md.dll )
      list( APPEND f90_system_lib libiomp5md.lib )
    endif( ENABLE_THREAD_PROFILE )
  else( ${MD_or_MT} MATCHES "MD" )
    if( ENABLE_THREAD_PROFILE )
      list( APPEND f90_system_lib libiompprof5mt.lib )
    else( ENABLE_THREAD_PROFILE )
      list( APPEND f90_system_lib libiomp5mt.lib )
    endif( ENABLE_THREAD_PROFILE )
  endif( ${MD_or_MT} MATCHES "MD" )
else( ${CMAKE_Fortran_COMPILER_VERSION_MAJOR} GREATER 10 )
  if( ${MD_or_MT} MATCHES "MD" )
    if( ENABLE_THREAD_PROFILE )
      list( APPEND f90_system_dll libguide40_stats.dll )
      list( APPEND f90_system_lib libguide40_stats.lib )
    else( ENABLE_THREAD_PROFILE )
      list( APPEND f90_system_dll libguide40.dll )
      list( APPEND f90_system_lib libguide40.lib )
    endif( ENABLE_THREAD_PROFILE )
  else( ${MD_or_MT} MATCHES "MD" )
    if( ENABLE_THREAD_PROFILE )
      list( APPEND f90_system_lib libguide_stats.lib )
    else( ENABLE_THREAD_PROFILE )
      list( APPEND f90_system_lib libguide.lib )
    endif( ENABLE_THREAD_PROFILE )
  endif( ${MD_or_MT} MATCHES "MD" )
endif( ${CMAKE_Fortran_COMPILER_VERSION_MAJOR} GREATER 10 )
  
  
   
# IVF dll directory
get_filename_component( CMAKE_Fortran_BIN_DIR ${CMAKE_Fortran_COMPILER} PATH )
# Useful for install commands (creating a portable developer version)
set( CMAKE_Fortran_redist_dll "" CACHE INTERNAL "Fortran redistributable system libraries that are needed by the applications built with Intel Visual Fortran." )
set( CMAKE_Fortran_debug_dll "" CACHE INTERNAL "Fortran system libraries that are needed by the applications built with Intel Visual Fortran (only optimized versions are redistributable.)" )
# Generate a list of F90 compiler libraries
foreach( lib ${f90_system_dll} )
  get_filename_component( libwe ${lib} NAME_WE )
  # optimized library
  find_file( CMAKE_Fortran_${libwe}_RELEASE
    NAMES ${libwe}.dll
    PATHS "${CMAKE_Fortran_BIN_DIR}"
    )
  mark_as_advanced( CMAKE_Fortran_${libwe}_RELEASE )
  set( CMAKE_Fortran_${libwe}_RELEASE "${CMAKE_Fortran_${libwe}_RELEASE}" 
       CACHE INTERNAL "Location of ${libwe}.dll for Intel Fortran." )
  # debug library
  find_file( CMAKE_Fortran_${libwe}_DEBUG
    NAMES 
       ${libwe}d.dll # debug
       ${libwe}.dll  # backup is opt version.
    PATHS "${CMAKE_Fortran_BIN_DIR}"
    )
  mark_as_advanced( CMAKE_Fortran_${libwe}_DEBUG )
  set( CMAKE_Fortran_${libwe}_DEBUG "${CMAKE_Fortran_${libwe}_DEBUG}" 
       CACHE INTERNAL "Location of ${libwe}.dll for Intel Fortran." )
  list( APPEND CMAKE_Fortran_redist_dll ${CMAKE_Fortran_${libwe}_RELEASE} )
  list( APPEND CMAKE_Fortran_debug_dll ${CMAKE_Fortran_${libwe}_DEBUG} )
  # set( CMAKE_Fortran_${libwe}_LIBRARY
    # "optimized;${CMAKE_Fortran_${libwe}_RELEASE};debug;${CMAKE_Fortran_${libwe}_DEBUG}"
    # CACHE INTERNAL "Fortran system libraries that are needed by the applications built with Intel Visual Fortran (only optimized versions are redistributable.)" FORCE )
endforeach( lib ${f90_system_dll} )

# Static libraries from the /Lib directory (useful for target_link_library command.
set( CMAKE_Fortran_compiler_libs "" CACHE INTERNAL "Fortran system libraries that are needed by the applications built with Intel Visual Fortran (only optimized versions are redistributable.)" )
string( REGEX REPLACE "[Bb]in" "lib"
         CMAKE_Fortran_LIB_DIR ${CMAKE_Fortran_BIN_DIR} )
if( NOT EXISTS ${CMAKE_Fortran_LIB_DIR} )
  message( FATAL_ERROR "Could not determine CMAKE_Fortran_LIB_DIR = ${CMAKE_Fortran_LIB_DIR}" )
endif( NOT EXISTS ${CMAKE_Fortran_LIB_DIR} )
         
foreach( lib ${f90_system_lib} )
  get_filename_component( libwe ${lib} NAME_WE )
  # optimized library
  find_file( CMAKE_Fortran_${libwe}_lib_RELEASE
    NAMES ${libwe}.lib
    PATHS "${CMAKE_Fortran_LIB_DIR}"
    NO_DEFAULT_PATH
    )
  mark_as_advanced( CMAKE_Fortran_${libwe}_lib_RELEASE )
  set( CMAKE_Fortran_${libwe}_lib_RELEASE "${CMAKE_Fortran_${libwe}_lib_RELEASE}" 
       CACHE INTERNAL "Location of ${libwe}.lib for Intel Fortran." )
  # debug library
  find_file( CMAKE_Fortran_${libwe}_lib_DEBUG
    NAMES 
       ${libwe}d.lib # debug
       ${libwe}.lib  # backup is opt version.
    PATHS "${CMAKE_Fortran_LIB_DIR}"
    NO_DEFAULT_PATH
    )
  mark_as_advanced( CMAKE_Fortran_${libwe}_lib_DEBUG )
  set( CMAKE_Fortran_${libwe}_lib_DEBUG "${CMAKE_Fortran_${libwe}_lib_DEBUG}" 
       CACHE INTERNAL "Location of ${libwe}.lib for Intel Fortran." )
  set( CMAKE_Fortran_${libwe}_lib_LIBRARY
    optimized
    "${CMAKE_Fortran_${libwe}_lib_RELEASE}"
    debug
    "${CMAKE_Fortran_${libwe}_lib_DEBUG}"
    CACHE INTERNAL "Fortran static system libraries that are needed by the applications built with Intel Visual Fortran (only optimized versions are redistributable.)" FORCE )
  list( APPEND CMAKE_Fortran_compiler_libs ${CMAKE_Fortran_${libwe}_lib_LIBRARY} )    
endforeach( lib ${f90_redist_lib} )
  
 
#------------------------------------------------------------------------------#
# End config/windows-ifort.cmake
#------------------------------------------------------------------------------#
