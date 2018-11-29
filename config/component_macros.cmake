#-----------------------------*-cmake-*----------------------------------------#
# file   config/component_macros.cmake
# author Kelly G. Thompson, kgt@lanl.gov
# date   2010 Dec 1
# brief  Provide extra macros to simplify CMakeLists.txt for component
#        directories.
# note   Copyright (C) 2016-2018 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------
# replacement for built in command 'add_executable'
#
# Purpose 1: In addition to adding an executable built from $sources, set
# Draco-specific properties for the library.  This macro reduces ~20
# lines of code repeated code down to 1-2.
#
# Purpose 2: Encapsulate library and vendor library dependencies per
# package.
#
# Purpose 3: Use information from 1 and 2 above to generate exported
# targets.
#
# Usage:
#
# add_component_executable(
#   TARGET       "target name"
#   EXE_NAME     "output executable name"
#   TARGET_DEPS  "dep1;dep2;..."
#   PREFIX       "ClubIMC"
#   SOURCES      "file1.cc;file2.cc;..."
#   HEADERS      "file1.hh;file2.hh;..."
#   VENDOR_LIST  "MPI;GSL"
#   VENDOR_LIBS  "${MPI_CXX_LIBRARIES};${GSL_LIBRARIES}"
#   VENDOR_INCLUDE_DIRS "${MPI_CXX_INCLUDE_DIR};${GSL_INCLUDE_DIR}"
#   FOLDER       "myfolder"
#   PROJECT_LABEL "myproject42"
#   NOEXPORT        - do not export target or dependencies to draco-config.cmake
#   NOCOMMANDWINDOW - On win32, do not create a command window (qt)
#   )
#
# Example:
#
# add_component_executable(
#   TARGET       Exe_draco_info
#   EXE_NAME     draco_info
#   TARGET_DEPS  Lib_diagnostics
#   PREFIX       Draco
#   SOURCES      "${PROJECT_SOURCE_DIR}/draco_info_main.cc"
#   FOLDER       diagnostics
#   )
#------------------------------------------------------------------------------
macro( add_component_executable )

  # These become variables of the form ${ace_NAME}, etc.
  cmake_parse_arguments(
    ace
    "NOEXPORT;NOCOMMANDWINDOW"
    "PREFIX;TARGET;EXE_NAME;LINK_LANGUAGE;FOLDER;PROJECT_LABEL"
    "HEADERS;SOURCES;TARGET_DEPS;VENDOR_LIST;VENDOR_LIBS;VENDOR_INCLUDE_DIRS"
    ${ARGV}
    )

  # Prefix for export
  if( NOT DEFINED ace_PREFIX AND NOT DEFINED ace_NOEXPORT)
    message( FATAL_ERROR
      "add_component_executable requires a PREFIX value to allow EXPORT of this target
or the target must be labeled NOEXPORT.")
  endif()

  # Default link language is C++
  if( "${ace_LINK_LANGUAGE}x" STREQUAL "x" )
    set( ace_LINK_LANGUAGE CXX )
  endif()

  #
  # Add headers to Visual Studio or Xcode solutions
  #
  if( ace_HEADERS )
    if( MSVC_IDE OR ${CMAKE_GENERATOR} MATCHES Xcode )
      list( APPEND ace_SOURCES ${ace_HEADERS} )
    endif()
  endif()

  if( NOT DEFINED ace_EXE_NAME )
    string( REPLACE "Exe_" "" ace_EXE_NAME ${ace_TARGET} )
  endif()

  #
  # Create the library and set the properties
  #

  # Set the component name: If registered from a test directory, extract
  # the parent's name.
  get_filename_component( ldir ${CMAKE_CURRENT_SOURCE_DIR} NAME )
  if( ${ldir} STREQUAL "test")
    get_filename_component( comp_target ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY )
    get_filename_component( comp_target ${comp_target} NAME )
    set( comp_target ${comp_target}_test )
  elseif( ${ldir} STREQUAL "bin" )
    get_filename_component( comp_target ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY )
    get_filename_component( comp_target ${comp_target} NAME )
  else()
    get_filename_component( comp_target ${CMAKE_CURRENT_SOURCE_DIR} NAME )
  endif()
  # Make the name safe: replace + with x
  string( REGEX REPLACE "[+]" "x" comp_target ${comp_target} )
  # Set the folder name:
  if( NOT DEFINED ace_FOLDER )
    set( ace_FOLDER ${comp_target} )
  endif()

  if( WIN32 AND ace_NOCOMMANDWINDOW )
    # The Win32 option prevents the command console from activating while the GUI is running.
    add_executable( ${ace_TARGET} WIN32 ${ace_SOURCES} )
  else()
    add_executable( ${ace_TARGET} ${ace_SOURCES} )
  endif()

  # Some properties are set at a global scope in compilerEnv.cmake:
  # - C_STANDARD, C_EXTENSIONS, CXX_STANDARD, CXX_EXTENSIONS,
  #   CXX_STANDARD_REQUIRED, and POSITION_INDEPENDENT_CODE
  set_target_properties( ${ace_TARGET} PROPERTIES
    OUTPUT_NAME ${ace_EXE_NAME}
    FOLDER      ${ace_FOLDER}
    INTERPROCEDURAL_OPTIMIZATION_RELEASE;${USE_IPO}
#    ENABLE_EXPORTS TRUE # See cmake policy cmp0065
    COMPILE_DEFINITIONS "PROJECT_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}\";PROJECT_BINARY_DIR=\"${PROJECT_BINARY_DIR}\""
    )
  if( DEFINED ace_PROJECT_LABEL )
    set_target_properties( ${ace_TARGET} PROPERTIES PROJECT_LABEL ${ace_PROJECT_LABEL} )
  endif()

  #
  # Generate properties related to library dependencies
  #
  if( DEFINED ace_TARGET_DEPS )
    target_link_libraries( ${ace_TARGET} ${ace_TARGET_DEPS} )
  endif()
  if( DEFINED ace_VENDOR_LIBS )
    target_link_libraries( ${ace_TARGET} ${ace_VENDOR_LIBS} )
  endif()
  if( DEFINED ace_VENDOR_INCLUDE_DIRS )
    include_directories( ${ace_VENDOR_INCLUDE_DIRS} )
  endif()

  #
  # Register the library for exported library support
  #

  # Find target file name and location
  get_target_property( impname ${ace_TARGET} OUTPUT_NAME )

  # the above command returns the location in the build tree.  We
  # need to convert this to the install location.
  # get_filename_component( imploc ${imploc} NAME )
  if( ${DRACO_SHARED_LIBS} )
    set( imploc "${CMAKE_INSTALL_PREFIX}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}${impname}${CMAKE_SHARED_LIBRARY_SUFFIX}" )
  else()
    set( imploc "${CMAKE_INSTALL_PREFIX}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}${impname}${CMAKE_STATIC_LIBRARY_SUFFIX}" )
  endif()

  set( ilil "")
  if( "${ace_TARGET_DEPS}x" STREQUAL "x" AND  "${ace_VENDOR_LIBS}x" STREQUAL "x")
    # do nothing
  elseif( "${ace_TARGET_DEPS}x" STREQUAL "x" )
    set( ilil "${ace_VENDOR_LIBS}" )
  elseif( "${ace_VENDOR_LIBS}x" STREQUAL "x")
    set( ilil "${ace_TARGET_DEPS}" )
  else()
    set( ilil "${ace_TARGET_DEPS};${ace_VENDOR_LIBS}" )
  endif()

  # For non-test libraries, save properties to the project-config.cmake file
  if( "${ilil}x" STREQUAL "x" )
    set( ${ace_PREFIX}_EXPORT_TARGET_PROPERTIES
      "${${ace_PREFIX}_EXPORT_TARGET_PROPERTIES}
    set_target_properties(${ace_TARGET} PROPERTIES
      IMPORTED_LINK_INTERFACE_LANGUAGES \"${ace_LINK_LANGUAGE}\"
      INTERFACE_INCLUDE_DIRECTORIES     \"${CMAKE_INSTALL_PREFIX}/${DBSCFGDIR}include\" )
    ")
  else()
    set( ${ace_PREFIX}_EXPORT_TARGET_PROPERTIES
      "${${ace_PREFIX}_EXPORT_TARGET_PROPERTIES}
    set_target_properties(${ace_TARGET} PROPERTIES
      IMPORTED_LINK_INTERFACE_LANGUAGES \"${ace_LINK_LANGUAGE}\"
      INTERFACE_INCLUDE_DIRECTORIES     \"${CMAKE_INSTALL_PREFIX}/${DBSCFGDIR}include\" )
  ")
  endif()

  # Only publish information to draco-config.cmake for non-test
  # libraries.  Also, omit any libraries that are marked as NOEXPORT
  if( NOT ${ace_NOEXPORT} AND
      NOT "${CMAKE_CURRENT_SOURCE_DIR}" MATCHES "test" )

    list( APPEND ${ace_PREFIX}_EXECUTABLES ${ace_TARGET} )
    list( APPEND ${ace_PREFIX}_TPL_LIST ${ace_VENDOR_LIST} )
    list( APPEND ${ace_PREFIX}_TPL_INCLUDE_DIRS ${ace_VENDOR_INCLUDE_DIRS} )
    list( APPEND ${ace_PREFIX}_TPL_LIBRARIES ${ace_VENDOR_LIBS} )
    if( ${ace_PREFIX}_TPL_INCLUDE_DIRS )
      list( REMOVE_DUPLICATES ${ace_PREFIX}_TPL_INCLUDE_DIRS )
    endif()
    if( ${ace_PREFIX}_TPL_LIBRARIES )
      list( REMOVE_DUPLICATES ${ace_PREFIX}_TPL_LIBRARIES )
    endif()
    if( ${ace_PREFIX}_TPL_LIST )
      list( REMOVE_DUPLICATES ${ace_PREFIX}_TPL_LIST )
    endif()
    if( ${ace_PREFIX}_EXECUTABLES )
      list( REMOVE_DUPLICATES ${ace_PREFIX}_EXECUTABLES )
    endif()

    set( ${ace_PREFIX}_EXECUTABLES "${${ace_PREFIX}_EXECUTABLES}"  CACHE INTERNAL "List of component targets" FORCE)
    set( ${ace_PREFIX}_TPL_LIST "${${ace_PREFIX}_TPL_LIST}"  CACHE INTERNAL
      "List of third party libraries known by ${ace_PREFIX}" FORCE)
    set( ${ace_PREFIX}_TPL_INCLUDE_DIRS "${${ace_PREFIX}_TPL_INCLUDE_DIRS}"  CACHE
      INTERNAL "List of include paths used by ${ace_PREFIX} to find thrid party vendor header files."
      FORCE)
    set( ${ace_PREFIX}_TPL_LIBRARIES "${${ace_PREFIX}_TPL_LIBRARIES}"  CACHE INTERNAL
      "List of third party libraries used by ${ace_PREFIX}." FORCE)
    set( ${ace_PREFIX}_EXPORT_TARGET_PROPERTIES
      "${${ace_PREFIX}_EXPORT_TARGET_PROPERTIES}" PARENT_SCOPE)

  endif()

  # If Win32, copy dll files into binary directory.
  copy_dll_link_libraries_to_build_dir( ${ace_TARGET} )

endmacro( add_component_executable )

#------------------------------------------------------------------------------
# replacement for built in command 'add_library'
#
# Purpose 1: In addition to adding a library built from $sources, set
# Draco-specific properties for the library.  This macro reduces ~20
# lines of code down to 1-2.
#
# Purpose 2: Encapsulate library and vendor library dependencies per
# package.
#
# Purpose 3: Use information from 1 and 2 above to generate exported
# targets.
#
# Usage:
#
# add_component_library(
#   TARGET       "target name"
#   LIBRARY_NAME "output library name"
#   TARGET_DEPS  "dep1;dep2;..."
#   PREFIX       "ClubIMC"
#   SOURCES      "file1.cc;file2.cc;..."
#   HEADERS      "file1.hh;file2.hh;..."
#   LIBRARY_NAME_PREFIX "rtt_"
#   VENDOR_LIST  "MPI;GSL"
#   VENDOR_LIBS  "${MPI_CXX_LIBRARIES};${GSL_LIBRARIES}"
#   VENDOR_INCLUDE_DIRS "${MPI_CXX_INCLUDE_DIR};${GSL_INCLUDE_DIR}"
#   NOEXPORT
#   )
#
# Example:
#
# add_component_library(
#   TARGET       Lib_quadrature
#   LIBRARY_NAME quadrature
#   TARGET_DEPS  "Lib_parser;Lib_special_functions;Lib_mesh_element"
#   PREFIX       "Draco"
#   SOURCES      "${sources}"
#   )
#
# Note: you must use quotes around ${list_of_sources} to preserve the list.
#------------------------------------------------------------------------------
macro( add_component_library )
  # target_name outputname sources
  # optional argument: libraryPrefix

  # These become variables of the form ${acl_NAME}, etc.
  cmake_parse_arguments(
    acl
    "NOEXPORT"
    "PREFIX;TARGET;LIBRARY_NAME;LIBRARY_NAME_PREFIX;LINK_LANGUAGE"
    "HEADERS;SOURCES;TARGET_DEPS;VENDOR_LIST;VENDOR_LIBS;VENDOR_INCLUDE_DIRS"
    ${ARGV}
    )

  #
  # Defaults:
  #
  # Optional 3rd argument is the library prefix.  The default is "rtt_".
  if( NOT acl_LIBRARY_NAME_PREFIX )
    set( acl_LIBRARY_NAME_PREFIX "rtt_" )
  endif()
  # Default link language is C++
  if( NOT acl_LINK_LANGUAGE )
    set( acl_LINK_LANGUAGE CXX )
  endif()

  #
  # Add headers to Visual Studio or Xcode solutions
  #
  if( acl_HEADERS )
    if( MSVC_IDE OR ${CMAKE_GENERATOR} MATCHES Xcode )
      list( APPEND acl_SOURCES ${acl_HEADERS} )
    endif()
  endif()

  #
  # Create the library and set the properties
  #

  # This is a test library.  Find the component name
  string( REPLACE "_test" "" comp_target ${acl_TARGET} )
  # extract project name, minus leading "Lib_"
  string( REPLACE "Lib_" "" folder_name ${acl_TARGET} )

  add_library( ${acl_TARGET} ${DRACO_LIBRARY_TYPE} ${acl_SOURCES} )
  # Some properties are set at a global scope in compilerEnv.cmake:
  # - C_STANDARD, C_EXTENSIONS, CXX_STANDARD, CXX_EXTENSIONS,
  #   CXX_STANDARD_REQUIRED, and POSITION_INDEPENDENT_CODE
  set_target_properties( ${acl_TARGET} PROPERTIES
    # ${compdefs}
    # Use custom library naming
    OUTPUT_NAME ${acl_LIBRARY_NAME_PREFIX}${acl_LIBRARY_NAME}
    FOLDER      ${folder_name}
    INTERPROCEDURAL_OPTIMIZATION_RELEASE;${USE_IPO}
    WINDOWS_EXPORT_ALL_SYMBOLS ON
    )

  #
  # Generate properties related to library dependencies
  #
  if( NOT "${acl_TARGET_DEPS}x" STREQUAL "x" )
    target_link_libraries( ${acl_TARGET} ${acl_TARGET_DEPS} )
  endif()
  if( NOT "${acl_VENDOR_LIBS}x" STREQUAL "x" )
    target_link_libraries( ${acl_TARGET} ${acl_VENDOR_LIBS} )
  endif()
  if( NOT "${acl_VENDOR_INCLUDE_DIRS}x" STREQUAL "x" )
    include_directories( ${acl_VENDOR_INCLUDE_DIRS} )
  endif()

  #
  # Register the library for exported library support
  #

  # Defaults
  if( "${acl_PREFIX}x" STREQUAL "x" )
    set( acl_PREFIX "Draco" )
  endif()

  # Find target file name and location
  get_target_property( impname ${acl_TARGET} OUTPUT_NAME )

  # the above command returns the location in the build tree.  We need to
  # convert this to the install location.
  if( ${DRACO_SHARED_LIBS} )
    set( imploc "${CMAKE_INSTALL_PREFIX}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}${impname}${CMAKE_SHARED_LIBRARY_SUFFIX}" )
  else()
    set( imploc "${CMAKE_INSTALL_PREFIX}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}${impname}${CMAKE_STATIC_LIBRARY_SUFFIX}" )
  endif()

  set( ilil "")
  if( "${acl_TARGET_DEPS}x" STREQUAL "x" AND "${acl_VENDOR_LIBS}x" STREQUAL "x")
    # do nothing
  elseif( "${acl_TARGET_DEPS}x" STREQUAL "x" )
    set( ilil "${acl_VENDOR_LIBS}" )
  elseif( "${acl_VENDOR_LIBS}x" STREQUAL "x")
    set( ilil "${acl_TARGET_DEPS}" )
  else()
    set( ilil "${acl_TARGET_DEPS};${acl_VENDOR_LIBS}" )
  endif()

  # For non-test libraries, save properties to the project-config.cmake file
  if( "${ilil}x" STREQUAL "x" )
    set( ${acl_PREFIX}_EXPORT_TARGET_PROPERTIES
      "${${acl_PREFIX}_EXPORT_TARGET_PROPERTIES}
    set_target_properties(${acl_TARGET} PROPERTIES
      IMPORTED_LINK_INTERFACE_LANGUAGES \"${acl_LINK_LANGUAGE}\"
      INTERFACE_INCLUDE_DIRECTORIES     \"${CMAKE_INSTALL_PREFIX}/${DBSCFGDIR}include\" )
    ")
  else()
#      IMPORTED_LINK_INTERFACE_LIBRARIES \"${ilil}\"
    set( ${acl_PREFIX}_EXPORT_TARGET_PROPERTIES
      "${${acl_PREFIX}_EXPORT_TARGET_PROPERTIES}
    set_target_properties(${acl_TARGET} PROPERTIES
      IMPORTED_LINK_INTERFACE_LANGUAGES \"${acl_LINK_LANGUAGE}\"
      INTERFACE_INCLUDE_DIRECTORIES     \"${CMAKE_INSTALL_PREFIX}/${DBSCFGDIR}include\" )
  ")
  endif()

  # Only publish information to draco-config.cmake for non-test libraries.
  # Also, omit any libraries that are marked as NOEXPORT
  if( NOT ${acl_NOEXPORT} AND
      NOT "${acl_TARGET}" MATCHES "test" )

    list( APPEND ${acl_PREFIX}_LIBRARIES ${acl_TARGET} )
    string( REPLACE "Lib_" "" compname ${acl_TARGET} )
    list( APPEND ${acl_PREFIX}_PACKAGE_LIST ${compname} )

    list( APPEND ${acl_PREFIX}_TPL_LIST ${acl_VENDOR_LIST} )
    list( APPEND ${acl_PREFIX}_TPL_INCLUDE_DIRS ${acl_VENDOR_INCLUDE_DIRS} )
    list( APPEND ${acl_PREFIX}_TPL_LIBRARIES ${acl_VENDOR_LIBS} )
    if( ${acl_PREFIX}_TPL_INCLUDE_DIRS )
      list( REMOVE_DUPLICATES ${acl_PREFIX}_TPL_INCLUDE_DIRS )
    endif()
    if( ${acl_PREFIX}_TPL_LIBRARIES )
      list( REMOVE_DUPLICATES ${acl_PREFIX}_TPL_LIBRARIES )
    endif()
    if( ${acl_PREFIX}_TPL_LIST )
      list( REMOVE_DUPLICATES ${acl_PREFIX}_TPL_LIST )
    endif()

    set( ${acl_PREFIX}_LIBRARIES "${${acl_PREFIX}_LIBRARIES}"  CACHE INTERNAL "List of component targets" FORCE)
    set( ${acl_PREFIX}_PACKAGE_LIST "${${acl_PREFIX}_PACKAGE_LIST}"  CACHE INTERNAL
      "List of known ${acl_PREFIX} targets" FORCE)
    set( ${acl_PREFIX}_TPL_LIST "${${acl_PREFIX}_TPL_LIST}"  CACHE INTERNAL
      "List of third party libraries known by ${acl_PREFIX}" FORCE)
    set( ${acl_PREFIX}_TPL_INCLUDE_DIRS "${${acl_PREFIX}_TPL_INCLUDE_DIRS}"  CACHE
      INTERNAL "List of include paths used by ${acl_PREFIX} to find thrid party vendor header files."
      FORCE)
    set( ${acl_PREFIX}_TPL_LIBRARIES "${${acl_PREFIX}_TPL_LIBRARIES}"  CACHE INTERNAL
      "List of third party libraries used by ${acl_PREFIX}." FORCE)
    set( ${acl_PREFIX}_EXPORT_TARGET_PROPERTIES
      "${${acl_PREFIX}_EXPORT_TARGET_PROPERTIES}" PARENT_SCOPE)

  endif()

endmacro()

# ------------------------------------------------------------
# Register_scalar_test()
#
# 1. Special treatment for Roadrunner/ppe code (must ssh and then run)
# 2. Register the test
# 3. Register the pass/fail criteria.
# ------------------------------------------------------------
macro( register_scalar_test targetname runcmd command cmd_args )

  separate_arguments( cmdargs UNIX_COMMAND ${cmd_args} )

  set(lverbose OFF)
  if( lverbose)
    message("add_test( NAME ${targetname} COMMAND ${RUN_CMD} ${command} ${cmdargs} )")
  endif()
  add_test( NAME ${targetname} COMMAND ${RUN_CMD} ${command} ${cmdargs} )

  # Reserve enough threads for application unit tests. Normally we only need 1
  # core for each scalar test.
  set( num_procs 1 )

  # For application unit tests, a parallel job is forked that needs more cores.
  if( addscalartest_APPLICATION_UNIT_TEST )
    if( "${cmd_args}" MATCHES "--np" AND NOT "${cmd_args}" MATCHES "scalar")
      string( REGEX REPLACE "--np ([0-9]+)" "\\1" num_procs "${cmd_args}" )
      # the forked processes needs $num_proc threads.  add one for the master
      # thread, the original scalar process.
      math( EXPR num_procs  "${num_procs} + 1" )
    endif()
  endif()

  # set pass fail criteria, processors required, etc.
  set_tests_properties( ${targetname}
    PROPERTIES
    PASS_REGULAR_EXPRESSION "${addscalartest_PASS_REGEX}"
    FAIL_REGULAR_EXPRESSION "${addscalartest_FAIL_REGEX}"
    PROCESSORS              "${num_procs}"
    WORKING_DIRECTORY       "${PROJECT_BINARY_DIR}"
    )
  if( NOT "${addscalartest_RESOURCE_LOCK}none" STREQUAL "none" )
    set_tests_properties( ${targetname}
      PROPERTIES RESOURCE_LOCK "${addscalartest_RESOURCE_LOCK}" )
  endif()
  if( NOT "${addscalartest_RUN_AFTER}none" STREQUAL "none" )
    set_tests_properties( ${targetname}
      PROPERTIES DEPENDS "${addscalartest_RUN_AFTER}" )
  endif()

  # Labels
  # message("LABEL (ast) = ${addscalartest_LABEL}")
  if( NOT "${addscalartest_LABEL}x" STREQUAL "x" )
    set_tests_properties( ${targetname}
      PROPERTIES  LABELS "${addscalartest_LABEL}" )
  endif()
  unset( num_procs )
  unset( lverbose )
endmacro()

# ------------------------------------------------------------
# Register_parallel_test()
#
# 1. Register the test
# 2. Register the pass/fail criteria.
# ------------------------------------------------------------
macro( register_parallel_test targetname numPE command cmd_args )
  set( lverbose OFF )
  if( lverbose )
    message( "      Adding test: ${targetname}" )
  endif()
  unset( RUN_CMD )

  if( addparalleltest_MPI_PLUS_OMP )
    string( REPLACE " " ";" mpiexec_omp_postflags_list "${MPIEXEC_OMP_POSTFLAGS}" )
    add_test(
      NAME    ${targetname}
      COMMAND ${RUN_CMD} ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${numPE}
              ${mpiexec_omp_postflags_list}
              ${command}
              ${cmdarg}
              )
  else()
    add_test(
      NAME    ${targetname}
      COMMAND ${RUN_CMD} ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${numPE}
              ${MPIRUN_POSTFLAGS}
              ${command}
              ${cmdarg}
              )
  endif()
  set_tests_properties( ${targetname}
    PROPERTIES
    PASS_REGULAR_EXPRESSION "${addparalleltest_PASS_REGEX}"
    FAIL_REGULAR_EXPRESSION "${addparalleltest_FAIL_REGEX}"
    WORKING_DIRECTORY       "${PROJECT_BINARY_DIR}"
    )
  if( NOT "${addparalleltest_RESOURCE_LOCK}none" STREQUAL "none" )
    set_tests_properties( ${targetname}
      PROPERTIES RESOURCE_LOCK "${addparalleltest_RESOURCE_LOCK}" )
  endif()
  if( NOT "${addparalleltest_RUN_AFTER}none" STREQUAL "none" )
    set_tests_properties( ${targetname}
      PROPERTIES DEPENDS "${addparalleltest_RUN_AFTER}" )
  endif()

  if( addparalleltest_MPI_PLUS_OMP )

    if( DEFINED ENV{OMP_NUM_THREADS} )
      math( EXPR numthreads "${numPE} * $ENV{OMP_NUM_THREADS}" )
    else()
      math( EXPR numthreads "${numPE} * ${MPI_CORES_PER_CPU}" )
    endif()

    if( MPI_HYPERTHREADING )
      math( EXPR numthreads "2 * ${numthreads}" )
    endif()
    set_tests_properties( ${targetname}
      PROPERTIES
        PROCESSORS "${numthreads}"
        LABELS     "nomemcheck" )
    unset( numthreads )
    unset( nnodes )
    unset( nnodes_remainder )

  else()

    if( DEFINED addparalleltest_LABEL )
      set_tests_properties( ${targetname}
          PROPERTIES LABELS "${addparalleltest_LABEL}" )
    endif()
    set_tests_properties( ${targetname} PROPERTIES PROCESSORS "${numPE}" )

  endif()
  unset( lverbose )
endmacro()

#----------------------------------------------------------------------#
# Special post-build options for Win32 platforms
# ---------------------------------------------------------------------#
# copy_dll_link_libraries_to_build_dir( target )
#
# For Win32 with shared libraries, all dll dependencies must be located
# in the PATH or in the application directory.  This cmake function
# creates POST_BUILD rules for unit tests and applications to ensure
# that the most up-to-date versions of all dependencies are in the same
# directory as the application.
#
# Consider replacing this functionality with CMake's BundleUtilities.
#----------------------------------------------------------------------#
function( copy_dll_link_libraries_to_build_dir target )

  if( NOT WIN32 )
    # Win32 platforms require all dll libraries to be in the local directory
    # (or $PATH)
    return()
  endif()

  # Debug dependencies for a particular target (uncomment the next line and
  # provide the targetname): "Ut_${compname}_${testname}_exe"
  if( "Ut_rng_ut_gsl_exe_foo" STREQUAL ${target} )
     set(lverbose ON)
  endif()
  if( lverbose )
     include(print_target_properties)
  endif()

  # For Win32 with shared libraries, the package dll must be located in the test
  # directory.

  # Discover all library dependencies for this unit test.
  get_target_property( link_libs ${target} LINK_LIBRARIES )
  if( lverbose )
    message("\nDebugging dependencies for target ${target}")
    # "${compname}_${testname}")
    message("  Dependencies = ${link_libs}\n")
  endif()
  if( "${link_libs}" MATCHES NOTFOUND )
     return() # nothing to do
  endif()

  set( old_link_libs "" )
  # Walk through the library dependencies to build a list of all .dll
  # dependencies.
  while( NOT "${old_link_libs}" STREQUAL "${link_libs}" )
    if(lverbose)
       message("Found new libraries (old_link_libs != link_libs).  Restarting search loop...\n")
    endif()
    set( old_link_libs ${link_libs} )
    foreach( lib ${link_libs} )
      if( lverbose )
        # message("\n  examine dependencies for lib           = ${lib}\n")
        print_targets_properties("${lib}")
      endif()
      # $lib will either be a cmake target (e.g.: Lib_dsxx, Lib_c4) or an actual
      # path to a library (c:\lib\gsl.lib).
      if( NOT EXISTS ${lib} AND TARGET ${lib} )
          # Must be a CMake target... find it's dependencies...
          # The target may be
          # 1. A target defined within the current build system (e.g.: Lib_c4), or
          # 2. an 'imported' targets like GSL::gsl.
          get_target_property( isimp ${lib} IMPORTED )
          if(isimp)
            get_target_property( link_libs3 ${lib} INTERFACE_LINK_LIBRARIES)
          else()
            get_target_property( link_libs2 ${lib} LINK_LIBRARIES )
          endif()
          list( APPEND link_libs ${link_libs2} )
          list( APPEND link_libs ${link_libs3} )
      endif()
    endforeach()
    # Loop through all current dependencies, remove static libraries
    #(they do not need to be in the run directory).
    list( REMOVE_DUPLICATES link_libs )
    foreach( lib ${link_libs} )
      if( "${lib}" MATCHES "NOTFOUND" )
        # nothing to add so remove from list
        list( REMOVE_ITEM link_libs ${lib} )
        if( lverbose )
          message("lib = ${lib} is NOTFOUND --> remove it from the list")
        endif()
      elseif( "${lib}" MATCHES "[$]<")
        # We have a generator expression.  This routine does not support this, so drop it.
        list( REMOVE_ITEM link_libs ${lib} )
        if( lverbose )
          message("lib = ${lib} is a generator expression --> remove it from the list")
        endif()
      elseif( "${lib}" MATCHES ".[lL]ib$" )
        # We have a path to a static library. Static libraries do not
        # need to be copied.
        list( REMOVE_ITEM link_libs ${lib} )
        if( lverbose )
          message("lib = ${lib} is a static lib --> remove it from the list")
        endif()
        # However, if there is a corresponding dll, we should add it
        # to the list.
        string( REPLACE ".lib" ".dll" dll_lib ${lib} )
        if( ${dll_lib} MATCHES "[.]dll$" AND EXISTS ${dll_lib} )
          list( APPEND link_libs "${dll_lib}" )
          if( lverbose )
            message("lib = ${lib} has a dll --> replace it with ${dll_lib}")
          endif()
        endif()
      endif()
    endforeach()
    if( lverbose )
      message("Updated dependencies list: ${target} --> ${link_libs}\n")
    endif()

  endwhile()

  list( REMOVE_DUPLICATES link_libs )
  # if( ${compname} MATCHES Fortran )
  # message("   ${compname}_${testname} --> ${link_libs}")
  # endif()

  if( lverbose )
    message("Creating post build commands for target = ${target}")
  endif()

  # Add a post-build command to copy each dll into the test directory.
  foreach( lib ${link_libs} )
    # We do not need to the post_build copy command for Draco Lib_* files.
    # These should already be in the correct location.
    if( NOT ${lib} MATCHES "Lib_" )
      if( lverbose )
        message("  looking at ${lib}")
      endif()
      unset( target_loc )
      if( NOT TARGET ${lib})
        if (lverbose )
          message("  ${lib} is not a target. skip it.")
        endif()
        continue()
      endif()
      # TYPE = {STATIC_LIBRARY, MODULE_LIBRARY, SHARED_LIBRARY,
      #         INTERFACE_LIBRARY, EXECUTABLE}
      # We cannot query INTERFACE_LIBRARY targets.
      get_target_property(lib_type ${lib} TYPE )
      if( ${lib_type} STREQUAL "INTERFACE_LIBRARY" )
        if( lverbose )
          message("  I think ${lib} is an INTERFACE_LIBRARY. Skipping to next dependency.")
        endif()
        continue()
      endif()
      get_target_property(is_imported ${lib} IMPORTED )
      if( is_imported )
        get_target_property(target_loc ${lib} IMPORTED_LOCATION_RELEASE )
        if( ${target_loc} MATCHES "NOTFOUND" )
          get_target_property(target_loc ${lib} IMPORTED_LOCATION_DEBUG )
        endif()
        if( ${target_loc} MATCHES "NOTFOUND" )
          get_target_property(target_loc ${lib} IMPORTED_LOCATION )
        endif()
        if( ${target_loc} MATCHES "NOTFOUND" )
          # path not found, ignore.
          if (lverbose )
            message("  ${lib} does not have an IMPORTED_LOCATION value. skip it.")
          endif()
          continue()
        endif()
        get_target_property(target_gnutoms ${lib} GNUtoMS)
      elseif( EXISTS ${lib} )
        # If $lib is a full path to a library, add it to the list
        if (lverbose )
          message("  ${lib} is the full path to a library; adding to the list.")
        endif()
        set( target_loc ${lib} )
        set( target_gnutoms NOTFOUND )
      else()
        if (lverbose )
          message("  ${lib} is a target that points to $<TARGET:${lib}>.")
        endif()
        set( target_loc $<TARGET_FILE:${lib}> )
        get_target_property( target_gnutoms ${lib} GNUtoMS )
      endif()

      add_custom_command( TARGET ${target} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${target_loc}
                ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${CMAKE_CFG_INTDIR} )

      if( lverbose )
        message("  CMAKE_COMMAND -E copy_if_different ${target_loc} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${CMAKE_CFG_INTDIR}")
      endif()

    endif()
  endforeach()
  unset( lverbose )
endfunction( copy_dll_link_libraries_to_build_dir )

#----------------------------------------------------------------------#
# add_scalar_tests
#
# Given a list of sources, create unit test executables, one exe for
# each source file.  Register the test to be run by ctest.
#
# Usage:
#
# add_scalar_tests(
#    SOURCES "${test_sources}"
#    [ DEPS    "${library_dependencies}" ]
#    [ TEST_ARGS     "arg1;arg2" ]
#    [ PASS_REGEX    "regex" ]
#    [ FAIL_REGEX    "regex" ]
#    [ RESOURCE_LOCK "lockname" ]
#    [ RUN_AFTER     "test_name" ]
#    [ APPLICATION_UNIT_TEST ]
# )
#
# Options:
#   APPLICATION_UNIT_TEST - (CI/CT only) If present, do not run the
#        test under 'aprun'.  ApplicationUnitTest based tests must be
#        run this way.  Setting this option when DRACO_C4==SCALAR will
#        reset any value provided in TEST_ARGS to be "--np scalar".
#   LINK_WITH_FORTRAN - Tell the compiler to use the Fortran compiler
#        for the final link of the test.  This is needed for Intel and
#        PGI.
#
#----------------------------------------------------------------------#
macro( add_scalar_tests test_sources )

  # These become variables of the form ${addscalartests_SOURCES}, etc.
  cmake_parse_arguments(
    addscalartest
    "APPLICATION_UNIT_TEST;LINK_WITH_FORTRAN;NONE"
    "LABEL"
    "DEPS;FAIL_REGEX;PASS_REGEX;RESOURCE_LOCK;RUN_AFTER;SOURCES;TEST_ARGS"
    ${ARGV}
    )

  # Special Cases:
  # ------------------------------------------------------------
  # On some platforms (Trinity), even scalar tests must be run underneath
  # MPIEXEC_EXECUTABLE (aprun):
  separate_arguments(MPIEXEC_POSTFLAGS)
  if( "${MPIEXEC_EXECUTABLE}" MATCHES "srun" )
    set( RUN_CMD ${MPIEXEC_EXECUTABLE} ${MPIEXEC_POSTFLAGS} -n 1 )
  else()
    unset( RUN_CMD )
  endif()

  # Special cases for tests that use the ApplicationUnitTest
  # framework (see c4/ApplicationUnitTest.hh).
  if( addscalartest_APPLICATION_UNIT_TEST )
    # This is a special case for Cray environments. For application unit tests,
    # the main test runs on the 'login' node (1 rank only) and the real test is
    # run under 'aprun'.  So we do not prefix the test command with 'aprun'.
    if( "${MPIEXEC_EXECUTABLE}" MATCHES "aprun" )
      unset( RUN_CMD )
    endif()

    # If this is an ApplicationUnitTest based test then the TEST_ARGS will look
    # like "--np 1;--np 2;--np 4".  For the case where DRACO_C4 = SCALAR, we
    # will automatically demote these arguments to "--np scalar."
    if( "${DRACO_C4}" MATCHES "SCALAR" )
      set( addscalartest_TEST_ARGS "--np scalar" )
    endif()

  endif()

  # Sanity Checks
  # ------------------------------------------------------------
  if( "${addscalartest_SOURCES}none" STREQUAL "none" )
    message( FATAL_ERROR "You must provide the keyword SOURCES and a list of sources when using the add_scalar_tests macro.  Please see draco/config/component_macros.cmake::add_scalar_tests() for more information." )
  endif()

  # Pass/Fail criteria
  if( "${addscalartest_PASS_REGEX}none" STREQUAL "none" )
    set( addscalartest_PASS_REGEX ".*[Tt]est: PASSED" )
  endif()
  if( "${addscalartest_FAIL_REGEX}none" STREQUAL "none" )
    set( addscalartest_FAIL_REGEX ".*[Tt]est: FAILED" )
    list( APPEND addscalartest_FAIL_REGEX ".*ERROR:.*" )
    list( APPEND addscalartest_FAIL_REGEX "forrtl: error" )
  endif()

  # Format resource lock command
  if( NOT "${addscalartest_RESOURCE_LOCK}none" STREQUAL "none" )
    set( addscalartest_RESOURCE_LOCK
      "RESOURCE_LOCK ${addscalartest_RESOURCE_LOCK}")
  endif()

  # What is the component name (always use Lib_${compname} as a dependency).
  string( REPLACE "_test" "" compname ${PROJECT_NAME} )

  # Loop over each test source files:
  # 1. Compile the executable
  # 2. Register the unit test

  # Generate the executable
  # ------------------------------------------------------------
  foreach( file ${addscalartest_SOURCES} )

    get_filename_component( testname ${file} NAME_WE )
    add_executable( Ut_${compname}_${testname}_exe ${file} )
    # Some properties are set at a global scope in compilerEnv.cmake:
    # - C_STANDARD, C_EXTENSIONS, CXX_STANDARD, CXX_EXTENSIONS,
    #   CXX_STANDARD_REQUIRED, and POSITION_INDEPENDENT_CODE
    set_target_properties( Ut_${compname}_${testname}_exe
      PROPERTIES
      OUTPUT_NAME ${testname}
      VS_KEYWORD  ${testname}
      FOLDER      ${compname}_test
      INTERPROCEDURAL_OPTIMIZATION_RELEASE;${USE_IPO}
      COMPILE_DEFINITIONS "PROJECT_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}\";PROJECT_BINARY_DIR=\"${PROJECT_BINARY_DIR}\""
      )
    # Do we need to use the Fortran compiler as the linker?
    if( addscalartest_LINK_WITH_FORTRAN )
      set_target_properties( Ut_${compname}_${testname}_exe
        PROPERTIES LINKER_LANGUAGE Fortran )
    endif()
    target_link_libraries(
      Ut_${compname}_${testname}_exe
      ${test_lib_target_name}
      ${addscalartest_DEPS}
      )

    # Special post-build options for Win32 platforms
    # ------------------------------------------------------------
    copy_dll_link_libraries_to_build_dir( Ut_${compname}_${testname}_exe )

  endforeach()

  # Register the unit test
  # ------------------------------------------------------------
  foreach( file ${addscalartest_SOURCES} )
    get_filename_component( testname ${file} NAME_WE )

    if( "${addscalartest_TEST_ARGS}none" STREQUAL "none" )
      register_scalar_test( ${compname}_${testname}
        "${RUN_CMD}" $<TARGET_FILE:Ut_${compname}_${testname}_exe> "" )
    else()
      set( iarg "0" )
      foreach( cmdarg ${addscalartest_TEST_ARGS} )
        math( EXPR iarg "${iarg} + 1" )
        register_scalar_test( ${compname}_${testname}_arg${iarg}
          "${RUN_CMD}" $<TARGET_FILE:Ut_${compname}_${testname}_exe> "${cmdarg}" )
      endforeach()
    endif()
  endforeach()

endmacro(add_scalar_tests)

#----------------------------------------------------------------------#
# add_parallel_tests
#
# Given a list of sources, create unit test executables, one exe for
# each source file.  Register the test to be run by ctest.
#
# Usage:
#
# add_parallel_tests(
#    SOURCES "${test_sources}"
#    DEPS    "${library_dependencies}"
#    PE_LIST "1;2;4" )
#
# Optional parameters that require arguments.
#
#    SOURCES         - semi-colon delimited list of files.
#    PE_LIST         - semi-colon deliminte list of integers (number
#                      of MPI ranks).
#    DEPS            - CMake target dependencies.
#    TEST_ARGS       - Command line arguments to use when running the test.
#    PASS_REGEX      - This regex must exist in the output to produce
#                      a 'pass.'
#    FAIL_REGEX      - If this regex exists in the output, the test
#                      will 'fail.'
#    RESOURCE_LOCK   - Tests with this common string identifier will
#                      not be run concurrently.
#    RUN_AFTER       - The argument to this option is a test name that
#                      must complete before the current test will be
#                      allowed to run
#    MPIFLAGS
#    LABEL           - Label that can be used to select tests via
#                      ctest's -R or -E options.
#
# Optional parameters that require arguments.
#
#    MPI_PLUS_OMP    - This bool indicates that the test uses OpenMP
#                      for each MPI rank.
#    LINK_WITH_FORTRAN - Use the Fortran compiler to perform the final
#                      link of the unit test.
#----------------------------------------------------------------------#
macro( add_parallel_tests )

  cmake_parse_arguments(
    addparalleltest
    "MPI_PLUS_OMP;LINK_WITH_FORTRAN"
    "LABEL"
    "DEPS;FAIL_REGEX;MPIFLAGS;PASS_REGEX;PE_LIST;RESOURCE_LOCK;RUN_AFTER;SOURCES;TEST_ARGS"
    ${ARGV}
    )

  set(lverbose OFF)

  # Sanity Check
  if( "${addparalleltest_SOURCES}none" STREQUAL "none" )
    message( FATAL_ERROR "You must provide the keyword SOURCES and a list of sources when using the add_parallel_tests macro.  Please see draco/config/component_macros.cmake::add_parallel_tests() for more information." )
  endif()

  # Pass/Fail criteria
  if( "${addparalleltest_PASS_REGEX}none" STREQUAL "none" )
    set( addparalleltest_PASS_REGEX ".*[Tt]est: PASSED" )
  endif()
  if( "${addparalleltest_FAIL_REGEX}none" STREQUAL "none" )
    set( addparalleltest_FAIL_REGEX ".*[Tt]est: FAILED" )
    list( APPEND addparalleltest_FAIL_REGEX ".*ERROR:.*" )
    list( APPEND addparalleltest_FAIL_REGEX "forrtl: error" )
  endif()

  # Format resource lock command
  if( NOT "${addparalleltest_RESOURCE_LOCK}none" STREQUAL "none" )
    set( addparalleltest_RESOURCE_LOCK
      "RESOURCE_LOCK ${addparalleltest_RESOURCE_LOCK}")
  endif()

  # What is the component name? Use this to give a target name to the test.
  string( REPLACE "_test" "" compname ${PROJECT_NAME} )

  # Override MPI Flags upon user request
  if ( NOT DEFINED addparalleltest_MPIFLAGS )
    set( MPIRUN_POSTFLAGS ${MPIEXEC_POSTFLAGS} )
  else()
    set( MPIRUN_POSTFLAGS "${addparalleltest_MPIFLAGS}" )
  endif()
  separate_arguments( MPIRUN_POSTFLAGS )

  # Loop over each test source files:
  # 1. Compile the executable
  # 2. Link against dependencies (libraries)

  foreach( file ${addparalleltest_SOURCES} )
    get_filename_component( testname ${file} NAME_WE )
    if( lverbose )
      message( "   add_executable( Ut_${compname}_${testname}_exe ${file} )
   set_target_properties(
      Ut_${compname}_${testname}_exe
      PROPERTIES
      OUTPUT_NAME ${testname}
      VS_KEYWORD  ${testname}
      FOLDER      ${compname}_test
      INTERPROCEDURAL_OPTIMIZATION_RELEASE;${USE_IPO}
      COMPILE_DEFINITIONS \"PROJECT_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}\";PROJECT_BINARY_DIR=\"${PROJECT_BINARY_DIR}\"\"
      )")
    endif()
    add_executable( Ut_${compname}_${testname}_exe ${file} )
    # Some properties are set at a global scope in compilerEnv.cmake:
    # - C_STANDARD, C_EXTENSIONS, CXX_STANDARD, CXX_EXTENSIONS,
    #   CXX_STANDARD_REQUIRED, and POSITION_INDEPENDENT_CODE
    set_target_properties(
      Ut_${compname}_${testname}_exe
      PROPERTIES
      OUTPUT_NAME ${testname}
      VS_KEYWORD  ${testname}
      FOLDER      ${compname}_test
      INTERPROCEDURAL_OPTIMIZATION_RELEASE;${USE_IPO}
      COMPILE_DEFINITIONS "PROJECT_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}\";PROJECT_BINARY_DIR=\"${PROJECT_BINARY_DIR}\""
      )
    if( addparalleltest_MPI_PLUS_OMP )
      if( ${CMAKE_GENERATOR} MATCHES Xcode )
        set_target_properties( Ut_${compname}_${testname}_exe
          PROPERTIES XCODE_ATTRIBUTE_ENABLE_OPENMP_SUPPORT YES )
      endif()
    endif()
    # Do we need to use the Fortran compiler as the linker?
    if( addparalleltest_LINK_WITH_FORTRAN )
      set_target_properties( Ut_${compname}_${testname}_exe
        PROPERTIES LINKER_LANGUAGE Fortran )
    endif()

    if( lverbose )
      message( "    target_link_libraries(
      Ut_${compname}_${testname}_exe
      ${test_lib_target_name}
      ${addparalleltest_DEPS} )")
    endif()
    target_link_libraries(
      Ut_${compname}_${testname}_exe
      ${test_lib_target_name}
      ${addparalleltest_DEPS}
      )

    # Special post-build options for Win32 platforms
    # ------------------------------------------------------------
    copy_dll_link_libraries_to_build_dir( Ut_${compname}_${testname}_exe )

  endforeach()

  # 3. Register the unit test
  # 4. Register the pass/fail criteria.
  if( ${DRACO_C4} MATCHES "MPI" )
    foreach( file ${addparalleltest_SOURCES} )
      get_filename_component( testname ${file} NAME_WE )
      foreach( numPE ${addparalleltest_PE_LIST} )
        set( iarg 0 )
        if( "${addparalleltest_TEST_ARGS}none" STREQUAL "none" )
          register_parallel_test(
            ${compname}_${testname}_${numPE}
            ${numPE}
            $<TARGET_FILE:Ut_${compname}_${testname}_exe>
            "" )
        else()
          foreach( cmdarg ${addparalleltest_TEST_ARGS} )
            math( EXPR iarg "${iarg} + 1" )
            register_parallel_test(
              ${compname}_${testname}_${numPE}_arg${iarg}
              ${numPE}
              $<TARGET_FILE:Ut_${compname}_${testname}_exe>
              ${cmdarg} )
          endforeach()
        endif()
      endforeach()
    endforeach()
  else()
    # DRACO_C4=SCALAR Mode:
    foreach( file ${addparalleltest_SOURCES} )
      set( iarg "0" )
      get_filename_component( testname ${file} NAME_WE )

      set( addscalartest_PASS_REGEX "${addparalleltest_PASS_REGEX}" )
      set( addscalartest_FAIL_REGEX "${addparalleltest_FAIL_REGEX}" )
      set( addscalartest_RESOURCE_LOCK "${addparalleltest_RESOURCE_LOCK}" )
      set( addscalartest_RUN_AFTER "${addparalleltest_RUN_AFTER}" )

      if( "${addparalleltest_TEST_ARGS}none" STREQUAL "none" )
        if( lverbose )
          message("   register_scalar_test( ${compname}_${testname}
          \"${RUN_CMD}\" $<TARGET_FILE:Ut_${compname}_${testname}_exe> \"\" )")
        endif()
        register_scalar_test( ${compname}_${testname}
          "${RUN_CMD}" $<TARGET_FILE:Ut_${compname}_${testname}_exe> "" )
      else()

        foreach( cmdarg ${addparalleltest_TEST_ARGS} )
          math( EXPR iarg "${iarg} + 1" )
          if( lverbose )
            message("   register_scalar_test( ${compname}_${testname}_arg${iarg}
            \"${RUN_CMD}\" ${testname} \"${cmdarg}\" ) ")
            endif()
          register_scalar_test( ${compname}_${testname}_arg${iarg}
            "${RUN_CMD}" $<TARGET_FILE:Ut_${compname}_${testname}_exe> "${cmdarg}" )
        endforeach()

      endif()

    endforeach()
  endif()

  unset( lverbose )
endmacro()

#------------------------------------------------------------------------------#
# provide_aux_files
#
# Call this macro from a package CMakeLists.txt to instruct the build system
# that some files should be copied from the source directory into the build
# directory.
# ------------------------------------------------------------------------------#
macro( provide_aux_files )

  cmake_parse_arguments(
    auxfiles
    "NONE"
    "SRC_EXT;DEST_EXT;FOLDER"
    "FILES;TARGETS"
    ${ARGV}
    )

  unset(required_files)
  foreach( file ${auxfiles_FILES} )
    get_filename_component( srcfilenameonly ${file} NAME )
    if( auxfiles_SRC_EXT )
      if( auxfiles_DEST_EXT )
        # replace SRC_EXT with DEST_EXT
        string( REPLACE "${auxfiles_SRC_EXT}" "${auxfiles_DEST_EXT}"
          srcfilenameonly "${srcfilenameonly}" )
      else()
        # strip SRC_EXT
        string( REPLACE ${auxfiles_SRC_EXT} "" srcfilenameonly
          ${srcfilenameonly} )
      endif()
    else()
      if( auxfiles_DEST_EXT )
        # add DEST_EXT
        set( srcfilenameonly "${srcfilenameonly}${auxfiles_DEST_EXT}" )
      endif()
    endif()
    set( outfile ${PROJECT_BINARY_DIR}/${srcfilenameonly} )
    add_custom_command(
      OUTPUT  ${outfile}
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${file} ${outfile}
      DEPENDS ${file}
      COMMENT "Copying ${file} to ${outfile}"
      )
    list( APPEND required_files "${outfile}" )
  endforeach()
  string( REPLACE "_test" "" compname ${PROJECT_NAME} )

  # Extra logic if multiple calls from the same directory.
  if( DEFINED Ut_${compname}_install_inputs_iarg )
    math( EXPR Ut_${compname}_install_inputs_iarg
      "${Ut_${compname}_install_inputs_iarg} + 1" )
  else()
    set( Ut_${compname}_install_inputs_iarg "0" CACHE INTERNAL
      "counter for each provide_aux_files command. Used to create individual
targets for copying support files.")
  endif()
  add_custom_target(
    Ut_${compname}_install_inputs_${Ut_${compname}_install_inputs_iarg}
    ALL
    DEPENDS ${required_files};${auxfiles_TARGETS}
    )
  if( auxfiles_FOLDER )
    set( folder_name ${auxfiles_FOLDER} )
  else()
    set( folder_name ${compname}_test )
  endif()
  set_target_properties(
    Ut_${compname}_install_inputs_${Ut_${compname}_install_inputs_iarg}
    PROPERTIES FOLDER ${folder_name}
    )

endmacro()

#------------------------------------------------------------------------------#
# PROCESS_AUTODOC_PAGES - Run configure_file(...) for all .dcc.in files found in
# the autodoc directory.  Destination will be the autodoc directory in the
# component binary directory.  The CMakeLists.txt in the draco/autodoc directory
# knows how to find these files.
#
# This allows CMAKE variables to be inserted into the .dcc files (e.g.:
# @Draco_VERSION@)
#
# E.g.: process_autodoc_pages()
#------------------------------------------------------------------------------#
macro( process_autodoc_pages )
  file( GLOB autodoc_in autodoc/*.in )
  foreach( file ${autodoc_in} )
    get_filename_component( dest_file ${file} NAME_WE )
    configure_file( ${file} ${PROJECT_BINARY_DIR}/autodoc/${dest_file}.dcc
      @ONLY )
  endforeach()
  file( GLOB images_in autodoc/*.jpg autodoc/*.png autodoc/*.gif )
  list( LENGTH images_in num_images )
  if( ${num_images} GREATER 0 )
    list( APPEND DOXYGEN_IMAGE_PATH "${PROJECT_SOURCE_DIR}/autodoc" )
  endif()
  set( DOXYGEN_IMAGE_PATH "${DOXYGEN_IMAGE_PATH}" CACHE PATH
     "List of directories that contain images for doxygen pages." FORCE )
  unset( images_in )
  unset( num_images )
endmacro()

#------------------------------------------------------------------------------#
# ADD_DIR_IF_EXISTS - A helper macro used for including sub-project directories
# from src/CMakeLists.
#------------------------------------------------------------------------------#
macro( add_dir_if_exists package )
  if( EXISTS ${PROJECT_SOURCE_DIR}/${package} )
    message( "   ${package}" )
    add_subdirectory( ${package} )
  endif()
endmacro()

#------------------------------------------------------------------------------#
# End config/component_macros.cmake
#------------------------------------------------------------------------------#
