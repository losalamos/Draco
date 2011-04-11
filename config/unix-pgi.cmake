#-----------------------------*-cmake-*----------------------------------------#
# file   config/unix-pgi.cmake
# author Kelly Thompson 
# date   2010 Nov 1
# brief  Establish flags for Linux64 - Intel C++
# note   Copyright © 2010 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

#
# Sanity Checks
# 

if( BUILD_SHARED_LIBS )
  message( FATAL_ERROR "Feature not available - yell at KT." )
endif( BUILD_SHARED_LIBS )

# Cannot use strict ansi flags on RedStorm
option( ENABLE_STRICT_ANSI "Turn on strict ANSI compliance?" ON )
if( ${ENABLE_STRICT_ANSI} )
   set( STRICT_ANSI_FLAGS "-Xa -A --no_using_std" )
endif()


#
# C++ libraries required by Fortran linker
# 

# execute_process(
#    COMMAND ${CMAKE_CXX_COMPILER} -show
#    TIMEOUT 5
#    RESULT_VARIABLE tmp
#    OUTPUT_VARIABLE pgiCC_show_output
#    ERROR_VARIABLE err
#    )
# string( REPLACE "\n" ";" pgiCC_show_output ${pgiCC_show_output} )
# foreach( line ${pgiCC_show_output} )
#    if( "${line}" MATCHES "COMPLIBOBJ" )
#       string( REGEX REPLACE ".*=" "" pgi_libdir ${line} )
#    endif()
# endforeach()
# message( STATUS "PGI Library Dir = ${pgi_libdir}")

#
# config.h settings
#

execute_process(
  COMMAND ${CMAKE_C_COMPILER} -V
  OUTPUT_VARIABLE ABS_C_COMPILER_VER
  )
string( REGEX REPLACE "Copyright.*" " " 
  ABS_C_COMPILER_VER ${ABS_C_COMPILER_VER} )
string( STRIP ${ABS_C_COMPILER_VER} ABS_C_COMPILER_VER )

execute_process(
  COMMAND ${CMAKE_CXX_COMPILER} -V
  OUTPUT_VARIABLE ABS_CXX_COMPILER_VER
  )
string( REGEX REPLACE "Copyright.*" " " 
  ABS_CXX_COMPILER_VER ${ABS_CXX_COMPILER_VER} )
string( STRIP ${ABS_CXX_COMPILER_VER} ABS_CXX_COMPILER_VER )

#
# Compiler Flags
# 

# Flags from Draco autoconf build system:
# -Xa
# -A                       ansi
# --no_using_std           Enable (disable) implicit use of the std
#                          namespace when standard header files are
#                          included. 
# --diag_suppress 940      Suppress warning #940
# --diag_suppress 11           "      "     # 11
# -DNO_PGI_OFFSET
# -Kieee                   Perform floating-point operations in strict
#                          conformance with the IEEE 754 standard. 
# --no_implicit_include    Disable implicit inclusion of source files
#                          as a method of finding definitions of
#                          template entities to be instantiated. 
# -Mdaz                    Enable (disable) mode to treat denormalized
#                          floating point numbers as zero.  -Mdaz is
#                          default for -tp p7-64 targets; -Mnodaz is
#                          default otherwise. 
# -pgf90libs               Link-time option to add the pgf90 runtime
#                          libraries, allowing mixed-language programming. 


if( CMAKE_GENERATOR STREQUAL "Unix Makefiles" )
  set( CMAKE_C_FLAGS                "-Kieee -Mdaz -pgf90libs" )
  set( CMAKE_C_FLAGS_DEBUG          "-g -O0") # -DDEBUG") 
  set( CMAKE_C_FLAGS_RELEASE        "-O3 -DNDEBUG" )
  set( CMAKE_C_FLAGS_MINSIZEREL     "${CMAKE_C_FLAGS_RELEASE}" )
  set( CMAKE_C_FLAGS_RELWITHDEBINFO "-O3 -DNDEBUG" )

  set( CMAKE_CXX_FLAGS                "${CMAKE_C_FLAGS} ${STRICT_ANSI_FLAGS} --no_implicit_include --diag_suppress 940 --diag_suppress 11 --diag_suppress 450 -DNO_PGI_OFFSET" )
  set( CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_C_FLAGS_DEBUG}")
  set( CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_C_FLAGS_RELEASE}")
  set( CMAKE_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_RELEASE}")
  set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}" )
ENDIF()

string( TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UPPER )

# Fortran libraries needed by the C++ linker when linking against
# LAPACK, etc.
# set( extra_libs pgc )# pgf90 pgf90_rpm1 pgf902 pgftnrtl pghpf2
# unset( PGI_EXTRA_F90_LIBS )
# message( STATUS "PGI_EXTRA_F90_LIBS = " )
# foreach( library ${extra_libs} )
#    find_library( lib_pgf90_lib${library}
#       NAMES ${library}
#       PATHS ${pgi_libdir}
#       )
#    list( APPEND PGI_EXTRA_F90_LIBS ${lib_pgf90_lib${library}} )
#    message( STATUS "     ${lib_pgf90_lib${library}}" )
# endforeach()

#------------------------------------------------------------------------------#
# End config/unix-pgi.cmake
#------------------------------------------------------------------------------#
