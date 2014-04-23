#.rst:
# CMakeAddFortranSubdirectory
# ---------------------------
# 
# Use a version of gfortran that is not available from within the current 
# project.  For example, use MinGW gfortran from Visual Studio if a 
# Fortran compiler is not found, or use GNU gfortran from a XCode/clang 
# build project.
#
# The 'add_fortran_subdirectory' function adds a subdirectory to a
# project that contains a Fortran only sub-project.  The module will
# check the current compiler and see if it can support Fortran.  If no
# Fortran compiler is found and the compiler is MSVC or if the Generator 
# is XCode, then this module will try to find a gfortran compiler in local 
# environment (e.g.: MinGW gfortran).  It will then use an external project 
# to build with alternate (MinGW/Unix) tools.  It will also create imported 
# targets for the libraries created.  
#
# For visual studio, this will only work if the Fortran code is built into a 
# dll, so BUILD_SHARED_LIBS is turned on in the project. In addition the 
# CMAKE_GNUtoMS option is set to on, so that the MS .lib files are created.  
#
# Usage is as follows:
#
# ::
#
#   cmake_add_fortran_subdirectory(
#    <subdir>                # name of subdirectory
#    PROJECT <project_name>  # project name in subdir top CMakeLists.txt
#    ARCHIVE_DIR <dir>       # dir where project places .lib files
#    RUNTIME_DIR <dir>       # dir where project places .dll files
#    LIBRARIES <lib>...      # names of library targets to import
#    LINK_LIBRARIES          # link interface libraries for LIBRARIES
#     [LINK_LIBS <lib> <dep>...]...
#    CMAKE_COMMAND_LINE ...  # extra command line flags to pass to cmake
#    NO_EXTERNAL_INSTALL     # skip installation of external project
#    )
#
# Relative paths in ARCHIVE_DIR and RUNTIME_DIR are interpreted with
# respect to the build directory corresponding to the source directory
# in which the function is invoked.
#
# Limitations:
#
# NO_EXTERNAL_INSTALL is required for forward compatibility with a
# future version that supports installation of the external project
# binaries during "make install".

#=============================================================================
# Copyright 2011-2012 Kitware, Inc.
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

set(_MS_MINGW_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR})
include(CheckLanguage)
include(ExternalProject)
include(CMakeParseArguments)

###--------------------------------------------------------------------------------####
function(_setup_mingw_config_and_build source_dir build_dir)

  # If MINGW_GFORTRAN is specified on the cmake command line, but the full path is not
  # provided, try to find the path.
  if( MINGW_GFORTRAN AND NOT EXISTS ${MINGW_GFORTRAN} )
    find_program(tmp_gfortran NAMES ${MINGW_GFORTRAN} )
    if( tmp_gfortran )
       set( MINGW_GFORTRAN "${tmp_gfortran}" )
    endif()
  endif()
 
  # Look for a MinGW gfortran.
  find_program(MINGW_GFORTRAN 
    NAMES ${MINGW_GFORTRAN} gfortran
    PATHS
      c:/MinGW/bin
      "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\Windows\\CurrentVersion\\Uninstall\\MinGW;InstallLocation]/bin"
    )  
  
  if(NOT EXISTS ${MINGW_GFORTRAN})
    message(FATAL_ERROR
      "gfortran not found, please install MinGW with the gfortran option."
      "Or set the cache variable MINGW_GFORTRAN to the full path. "
      " This is required to build")
  endif()

  # Validate flavor/architecture of specified gfortran 
  if( MSVC )
      # MinGW gfortran under MSVS.
      if(CMAKE_SIZEOF_VOID_P EQUAL 8)
        set(_mingw_target "Target:.*64.*mingw")
      else()
        set(_mingw_target "Target:.*mingw32")
      endif()
  else() 
      # GNU gfortran under XCode.
      if(CMAKE_SIZEOF_VOID_P EQUAL 8)
        set(_mingw_target "Target:.*64-apple*")
      else()
        set(_mingw_target "Target:.*86-apple*")
      endif()  
  endif() # MSVC
  execute_process(COMMAND "${MINGW_GFORTRAN}" -v
    ERROR_VARIABLE out ERROR_STRIP_TRAILING_WHITESPACE)
  if(NOT "${out}" MATCHES "${_mingw_target}")
    string(REPLACE "\n" "\n  " out "  ${out}")
    message(FATAL_ERROR
      "MINGW_GFORTRAN is set to\n"
      "  ${MINGW_GFORTRAN}\n"
      "which is not a MinGW gfortran for this architecture.  "
      "The output from -v does not match \"${_mingw_target}\":\n"
      "${out}\n"
      "Set MINGW_GFORTRAN to a proper MinGW gfortran for this architecture."
      )
  endif()
      
  # Configure scripts to run MinGW tools with the proper PATH.
  get_filename_component(MINGW_PATH ${MINGW_GFORTRAN} PATH)
  file(TO_NATIVE_PATH "${MINGW_PATH}" MINGW_PATH)
  string(REPLACE "\\" "\\\\" MINGW_PATH "${MINGW_PATH}")
  # Generator type
  if( MSVC )
    set( GENERATOR_TYPE "-GMinGW Makefiles")
  else() # XCode
    set( GENERATOR_TYPE "-GUnix Makefiles")
  endif()
  configure_file(
    ${_MS_MINGW_SOURCE_DIR}/CMakeAddFortranSubdirectory/config_mingw.cmake.in
    ${build_dir}/config_mingw.cmake
    @ONLY)
  configure_file(
    ${_MS_MINGW_SOURCE_DIR}/CMakeAddFortranSubdirectory/build_mingw.cmake.in
    ${build_dir}/build_mingw.cmake
    @ONLY)
endfunction()

###--------------------------------------------------------------------------------####
function(_add_fortran_library_link_interface library depend_library)
  set_target_properties(${library} PROPERTIES
    IMPORTED_LINK_INTERFACE_LIBRARIES_NOCONFIG "${depend_library}")
endfunction()

###--------------------------------------------------------------------------------####
function(cmake_add_fortran_subdirectory subdir)
  # Parse arguments to function
  set(options NO_EXTERNAL_INSTALL)
  set(oneValueArgs PROJECT ARCHIVE_DIR RUNTIME_DIR)
  set(multiValueArgs LIBRARIES TARGET_NAMES LINK_LIBRARIES CMAKE_COMMAND_LINE)
  cmake_parse_arguments(ARGS "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  if(NOT ARGS_NO_EXTERNAL_INSTALL)
    message(FATAL_ERROR
      "Option NO_EXTERNAL_INSTALL is required (for forward compatibility) "
      "but was not given." )
  endif()

  # if we are not using MSVC without Fortran support then just use the usual 
  # add_subdirectory to build the Fortran library:
  check_language(Fortran)  
  if( CMAKE_Fortran_COMPILER )
    add_subdirectory(${subdir})
    return()
  endif()
  
  if( NOT (MSVC OR ${CMAKE_GENERATOR} MATCHES Xcode) )
     message( FATAL_ERROR "Add_fortran_subdirectory only tested for MSVC and XCode." )
  endif()

  # if we have MSVC without Intel Fortran then setup external projects to 
  # build with alternate Fortran

  set(source_dir   "${CMAKE_CURRENT_SOURCE_DIR}/${subdir}")
  set(project_name "${ARGS_PROJECT}")
  set(library_dir  "${ARGS_ARCHIVE_DIR}")
  set(binary_dir   "${ARGS_RUNTIME_DIR}")
  set(libraries     ${ARGS_LIBRARIES})
  set(target_names "${ARGS_TARGET_NAMES}")
  list(LENGTH libraries numlibs)
  list(LENGTH target_names numtgtnames)
  if( ${numtgtnames} STREQUAL 0 )
     set(target_names ${libraries})
     set( numtgtnames ${numlibs})
  endif()
  if( NOT ${numlibs} STREQUAL ${numtgtnames} )
     message(FATAL_ERROR "If TARGET_NAMES are provided, you must provide an "
     "equal number of entries for both TARGET_NAMES and LIBRARIES." )
  endif()
  # use the same directory that add_subdirectory would have used
  set(build_dir "${CMAKE_CURRENT_BINARY_DIR}/${subdir}")
  foreach(dir_var library_dir binary_dir)
    if(NOT IS_ABSOLUTE "${${dir_var}}")
      get_filename_component(${dir_var}
        "${CMAKE_CURRENT_BINARY_DIR}/${${dir_var}}" ABSOLUTE)
    endif()
  endforeach()
  # create build and configure wrapper scripts
  _setup_mingw_config_and_build("${source_dir}" "${build_dir}")
  # create the external project
  externalproject_add(${project_name}_build
    SOURCE_DIR ${source_dir}
    BINARY_DIR ${build_dir}
    CONFIGURE_COMMAND ${CMAKE_COMMAND}
    -P ${build_dir}/config_mingw.cmake
    BUILD_COMMAND ${CMAKE_COMMAND}
    -P ${build_dir}/build_mingw.cmake
    INSTALL_COMMAND ""
    )
  # make the external project always run make with each build
  externalproject_add_step(${project_name}_build forcebuild
    COMMAND ${CMAKE_COMMAND}
    -E remove
    ${CMAKE_CURRENT_BUILD_DIR}/${project_name}-prefix/src/${project_name}-stamp/${project_name}-build
    DEPENDEES configure
    DEPENDERS build
    ALWAYS 1
    )
  # create imported targets for all libraries
  set(idx 0)
  foreach(lib ${libraries})
    list(GET target_names idx tgt)
    add_library(${tgt} SHARED IMPORTED GLOBAL)
    set_property(TARGET ${tgt} APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
    set_target_properties(${tgt} PROPERTIES
      IMPORTED_IMPLIB_NOCONFIG   "${library_dir}/lib${lib}${CMAKE_STATIC_LIBRARY_SUFFIX}" #.LIB
      IMPORTED_LOCATION_NOCONFIG "${binary_dir}/lib${lib}${CMAKE_SHARED_LIBRARY_SUFFIX}"  #.DLL
      )
    add_dependencies(${tgt} ${project_name}_build)
  endforeach()

  # now setup link libraries for targets
  set(start FALSE)
  set(target)
  foreach(lib ${ARGS_LINK_LIBRARIES})
    if("${lib}" STREQUAL "LINK_LIBS")
      set(start TRUE)
    else()
      if(start)
        if(DEFINED target)
          # process current target and target_libs
          _add_fortran_library_link_interface(${target} "${target_libs}")
          # zero out target and target_libs
          set(target)
          set(target_libs)
        endif()
        # save the current target and set start to FALSE
        set(target ${lib})
        set(start FALSE)
      else()
        # append the lib to target_libs
        list(APPEND target_libs "${lib}")
      endif()
    endif()
  endforeach()
  # process anything that is left in target and target_libs
  if(DEFINED target)
    _add_fortran_library_link_interface(${target} "${target_libs}")
  endif()
endfunction()
