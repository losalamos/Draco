#-----------------------------*-cmake-*----------------------------------------#
# file   config/compiler_env.cmake
# author Kelly Thompson
# date   2008 May 30
# brief  Default CMake build parameters
# note   Copyright © 2010 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

# Global

# Use SSE{1,2,3} instruction set?
option( DRACO_ENABLE_SSE "Try to use SSE instructions?" ON )
mark_as_advanced( DRACO_ENABLE_SSE )

option( ENABLE_OPENMP "Link against OpenMP libraries?" OFF )
mark_as_advanced( ENABLE_OPENMP )

# Library type to build
# Linux: STATIC is a lib<XXX>.a
#        SHARED is a lib<XXX>.so (requires rpath or .so found in $LD_LIBRARY_PATH
# MSVC : STATIC is <XXX>.lib
#        SHARED is <XXX>.dll (requires dll to be in $PATH or in same directory as exe).
set( DRACO_LIBRARY_TYPE "SHARED" CACHE STRING 
	"Keyword for creating new libraries (STATIC or SHARED).")

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
  
  set( gen_comp_env_set 1 )
endmacro()

#------------------------------------------------------------------------------#
# Setup C++ Compiler
#------------------------------------------------------------------------------#
macro(dbsSetupCxx)
  
  if( NOT gen_comp_env_set STREQUAL 1 )
   dbsSetupCompilers()
  endif()
  
  if( ${CMAKE_CXX_COMPILER} MATCHES "cl" )
    include( windows-cl )
  elseif( ${CMAKE_CXX_COMPILER} MATCHES "ppu-g[+][+]" )
     include( unix-ppu )
  elseif( ${CMAKE_CXX_COMPILER} MATCHES "icpc" )
     include( unix-intel )
  elseif( ${CMAKE_CXX_COMPILER} MATCHES "pgCC" )
     include( unix-pgi )
  elseif( ${CMAKE_CXX_COMPILER} MATCHES "xt-asyncpe" ) # Ceilo (catamount/pgi)
     include( unix-pgi )
  elseif( ${CMAKE_CXX_COMPILER} MATCHES "[cg][+]+" )
    include( unix-g++ )
  else( ${CMAKE_CXX_COMPILER} MATCHES "cl" )
    message( FATAL_ERROR "Build system does not support CXX=${CMAKE_CXX_COMPILER}" )
  endif( ${CMAKE_CXX_COMPILER} MATCHES "cl" )
 
endmacro()

#------------------------------------------------------------------------------#
# Setup Fortran Compiler
#
# Use:
#    include( compilerEnv )
#    dbsSetupF90( [QUIET] )
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
#    ABS_FLOAT_PRECISION     - string (config.h)
#    PRECISION_DOUBLE | PRECISION_SINGLE - bool
# 
#------------------------------------------------------------------------------#
macro(dbsSetupFortran)
  
  # if( NOT gen_comp_env_set STREQUAL 1 )
  #  dbsSetupCompilers()
  # endif()
  
  if( ${CMAKE_Fortran_COMPILER} MATCHES "gfortran" )
    include( unix-gfortran )
#  elseif( ${CMAKE_Fortran_COMPILER} MATCHES "ifort" )
#    include( unix-ifort )
#  elseif( ${CMAKE_CXX_COMPILER} MATCHES "pgf90" )
#     include( unix-pgf90 )
  else()
    message( FATAL_ERROR "Build system does not support F90=${CMAKE_Fortran_COMPILER}" )
  endif()
 
endmacro()


# macro(setup_f90compiler)

  # if( ${ARGV} MATCHES QUIET )
    # set( QUIET "QUIET" )
  # endif( ${ARGV} MATCHES QUIET )

  # if( NOT gen_comp_env_set STREQUAL 1 )
   # setup_generic_compiler_env()
  # endif( NOT gen_comp_env_set STREQUAL 1 )

  # Tell CMake to look for a Fortran compiler and setup default flags.
#  enable_language(Fortran)

  # Do we use doubles or floats to represent Fortran floating point numbers?
  # option( ENABLE_SINGLE_PRECISION
    # "Compile the library with single precision floating point" ON )
  # if( ENABLE_SINGLE_PRECISION )
    # set( ABS_FLOAT_PRECISION single )
    # set( PRECISION_SINGLE ON )
  # else( ENABLE_SINGLE_PRECISION )
    # set( ABS_FLOAT_PRECISION double )
    # set( PRECISION_DOUBLE ON )
  # endif( ENABLE_SINGLE_PRECISION )

  # Setup compiler specific flags.
  # if( WIN32 )
    # if( ${CMAKE_Fortran_COMPILER} MATCHES "ifort" )
      # set_mt_or_md()
      # include( windows-ifort )
    # elseif( ${CMAKE_Fortran_COMPILER} MATCHES "g95" )
      # include( windows-g95 )
    # else( ${CMAKE_Fortran_COMPILER} MATCHES "ifort" )
      # message( FATAL_ERROR "Build system does not support FC=${CMAKE_Fortran_COMPILER}" )
    # endif( ${CMAKE_Fortran_COMPILER} MATCHES "ifort" )
  # elseif( UNIX )
    # if( ${CMAKE_Fortran_COMPILER} MATCHES "g95" )
      # include( unix-g95 )
    # elseif( ${CMAKE_Fortran_COMPILER} MATCHES "pathf9[05]" )
      # message( FATAL_ERROR "Build system does not support FC=${CMAKE_Fortran_COMPILER}" )
      # include( unix-pathf90 )
    # elseif( ${CMAKE_Fortran_COMPILER} MATCHES "gfortran" )
      # include( unix-gfortran )
    # elseif(  ${CMAKE_Fortran_COMPILER} MATCHES "ifort" )
      # include( unix-ifort )
    # else( ${CMAKE_Fortran_COMPILER} MATCHES "g95" )
      # message( FATAL_ERROR "Build system does not support FC=${CMAKE_Fortran_COMPILER}" )
    # endif( ${CMAKE_Fortran_COMPILER} MATCHES "g95" )
  # else( WIN32 )
    # message( FATAL_ERROR "Unsupported platform (not WIN32 and not UNIX )." )
  # endif( WIN32 )  

  # Determine how this Fortran compiler mangles names in libraries.
  #include( determineFortranNameFormats )
  #determineFunctionNameMangling( ${QUIET} )

  # Determine the naming scheme this Fortran compiler uses for .mod files.
  # determineModuleName( ${QUIET} )
  
  # if( NOT QUIET )
    # message("ENABLE_SINGLE_PRECISION : ${ENABLE_SINGLE_PRECISION}")
    # message("ABS_FLOAT_PRECISION     : ${ABS_FLOAT_PRECISION}")
    # message("PRECISION_SINGLE        : ${PRECISION_SINGLE}")
    # message("PRECISION_DOUBLE        : ${PRECISION_DOUBLE}")
    # message("ENABLE_OPENMP           : ${ENABLE_OMPENMP}")
    # message("ENABLE_THREAD_PROFILE   : ${ENABLE_THREAD_PROFILE}")
  # endif( NOT QUIET )
  
# endmacro(setup_f90compiler)


#======================================================================
# Predict the resulting module filename when the F90 source filename 
# with module ${modname} is known.
#
# Output: ${modfile} 
#======================================================================
# macro( mangle_module_name sourcefilename modulename )

  # set( sourcefilename ${ARGV0} )
  # set( modulename ${ARGV1} )
  
  # if( NOT CMAKE_Fortran_have_mod_scheme )
    # message( FATAL_ERROR "Unable to determine naming scheme for Fortran module files." )
  # endif( NOT CMAKE_Fortran_have_mod_scheme )
  
  # if( ${CMAKE_Fortran_mod_scheme} MATCHES MODULENAME.mod ) # Unix-pathscale
    ##moduleSRCFILE.f90,"module moduleNAME" --> MODULENAME.mod
    # string( TOUPPER ${modulename} tmp )
    # set( modfile "${tmp}.mod" )
  # elseif( ${CMAKE_Fortran_mod_scheme} MATCHES modulename.mod ) # Unix-g95
    ##moduleSRCFILE.f90,"module moduleNAME" --> modulename.mod
    # string( TOLOWER ${modulename} tmp )
    # set( modfile "${tmp}.mod" )    
  # elseif( ${CMAKE_Fortran_mod_scheme} MATCHES moduleNAME.mod ) # Windows-ifort
    ##moduleSRCFILE.f90,"module moduleNAME" --> moduleNAME.mod
    # set( modfile "${modulename}.mod" )
  # else( ${CMAKE_Fortran_mod_scheme} MATCHES MODULENAME.mod )
    # message( FATAL_ERROR "CMake script needs to be extended to support CMAKE_Fortran_mod_scheme = ${CMAKE_Fortran_mod_scheme}." )
  # endif( ${CMAKE_Fortran_mod_scheme} MATCHES MODULENAME.mod )

  ##return ${modfile}
  
# endmacro( mangle_module_name sourcefilename modulename )

#======================================================================
# Specify the name of the created modfile:
#
# Output: ${modfile_name} 
# Input : fsrc - An existing Fortran file.
#======================================================================
# macro( fortran_mod_name fsrc )

  ##intput variable
  # set( fsrc ${ARGV1} )

  ##output variable
  # set( modfile_name NOTFOUND )

  ##strip path from source filename.
  # get_filename_component( sourcefilename ${fsrc} NAME )

  ##name of module (not necessariily the filename!
  # if( NOT EXISTS "${fsrc}" )
    # message( FATAL_ERROR "Could not find source file ${fsrc} when trying to determine resulting module filename." )
  # endif( NOT EXISTS "${fsrc}" )

  ##To determine if this Fortran source will module file, we read the
  ##file and look for the key words "module ..."
  # file( READ ${fsrc} file_contents ) 
  ##Convert file contents string into individual lines.
  # string( REGEX REPLACE "\n" ";" lines "${file_contents}" )
  ##Look at each line looking for "module ..."
  # foreach( line ${lines} )
    ##"module" can be capital, lc or mixed.
    ##"module" must be proceded only by whitespace.
    ##name of module must consist of letters, numbers or underscores
    # if( ${line} MATCHES "^[ \t]*[Mm][Oo][Dd][Uu][Ll][Ee][ \t]+[A-z0-9_]+" )

      ##Extract the name of the module only.
      # string( REGEX REPLACE "^[ \t]*[Mm][Oo][Dd][Uu][Ll][Ee][ \t]+([A-z0-9_]+)" "\\1" modname "${line}" ) 
      # mangle_module_name( ${sourcefilename} ${modname} ) # returns ${modfile}
      # if( ${modfile_name} MATCHES NOTFOUND )
        # set( modfile_name ${PROJECT_BINARY_DIR}/${modfile} )
      # elseif( ${modfile_name} MATCHES NOTFOUND )
        # message( FATAL_ERROR "
# Build system assumes only one module produced per source file.  However, we
# found at least two module definitions:
# source file  = ${fcrc}
# first module = ${modfile_name}
# 2nd module   = ${PROJECT_BINARY_DIR}/${modfile}
# ")
      # endif( ${modfile_name} MATCHES NOTFOUND )

     # endif( ${line} MATCHES "^[ \t]*[Mm][Oo][Dd][Uu][Ll][Ee][ \t]+[A-z0-9_]+" )

  # endforeach( line ${lines} )

  ##set parent scope for $modfile_name
  # set( modfile_name "${modfile_name}" PARENT_SCOPE )

#endmacro( fortran_mod_name )

#======================================================================
# Create a list of FQ modfiles
# 
# modfiles will only contain a list of files that contain Fortran
# modules 
#
# Input : A list of 1 or more fortran source file names.
# Output: ${modfiles} 
#======================================================================
# macro( create_modfiles_fpp_src )

  # set( src_F ${ARGV} )
  # set( modfiles )
  
  # foreach( src_F ${F_src} )
 
    # set( src "${PROJECT_SOURCE_DIR}/${src_F}" )
    # fortran_mod_name( ${src} ) # returns ${modfile_name}
    # if( NOT ${modfile_name} MATCHES NOTFOUND )
      # list( APPEND modfiles "${modfile_name}" )
    # endif( NOT ${modfile_name} MATCHES NOTFOUND )

  # endforeach( src_F ${F_src} )
  
  ##extra clean items
  # set_directory_properties(
    # PROPERTIES 
    # ADDITIONAL_MAKE_CLEAN_FILES "${modfiles}"
    # )

# endmacro( create_modfiles_fpp_src )

#------------------------------------------------------------------------------#
# End config/compiler_env.cmake
#------------------------------------------------------------------------------#
