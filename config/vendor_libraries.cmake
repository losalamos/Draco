#-----------------------------*-cmake-*----------------------------------------#
# file   config/vendor_libraries.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 June 6
# brief  Look for any libraries which are required at the top level.
# note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

include( FeatureSummary )
include( setupMPI ) # defines the macros setupMPILibrariesUnix|Windows

#------------------------------------------------------------------------------
# Helper macros for LAPACK/Unix
#
# This module sets the following variables:
# lapack_FOUND - set to true if a library implementing the LAPACK
#         interface is found
# lapack_VERSION - '3.4.1'
# provides targets: lapack, blas
#------------------------------------------------------------------------------
macro( setupLAPACKLibrariesUnix )

  # There are several flavors of LAPACK.
  # 1. look for netlib-lapack
  # 2. look for MKL (Intel)
  # 3. look for OpenBLAS.

  message( STATUS "Looking for lapack (netlib)...")
  set( lapack_FOUND FALSE )

  # Use LAPACK_LIB_DIR, if the user set it, to help find LAPACK.
  # This first try will also look for BLAS/LAPACK at CMAKE_PREFIX_PATH.
  if( EXISTS ${LAPACK_LIB_DIR}/cmake )
    file( GLOB lapack_cmake_prefix_path
      LIST_DIRECTORIES true
      ${LAPACK_LIB_DIR}/cmake/lapack-* )
    list( APPEND CMAKE_PREFIX_PATH ${lapack_cmake_prefix_path} )
  endif()
  find_package( lapack CONFIG QUIET )

  if( lapack_FOUND )
    set( lapack_flavor "netlib")
    foreach( config NOCONFIG DEBUG RELEASE RELWITHDEBINFO )
      get_target_property(tmp lapack IMPORTED_LOCATION_${config} )
      if( EXISTS ${tmp} )
        set( lapack_FOUND TRUE )
      endif()
    endforeach()
    message( STATUS "Looking for lapack (netlib)....found ${LAPACK_LIB_DIR}")
    set( lapack_FOUND ${lapack_FOUND} CACHE BOOL "Did we find LAPACK." FORCE )

    # The above might define blas, or it might not. Double check:
    if( NOT TARGET blas )
      find_package( BLAS )
      if( BLAS_FOUND )
        add_library( blas STATIC IMPORTED)
        set_target_properties( blas PROPERTIES
          IMPORTED_LOCATION                 "${BLAS_LIBRARIES}"
          IMPORTED_LINK_INTERFACE_LANGUAGES "Fortran")
      else()
        message( FATAL_ERROR "Looking for lapack (netlib)....blas not found")
      endif()
    endif()

  else()
    message( STATUS "Looking for lapack (netlib)....not found")
  endif()

  mark_as_advanced( lapack_DIR lapack_FOUND )

  # Above we tried to find lapack-config.cmake at $LAPACK_LIB_DIR/cmake/lapack.
  # This is a draco supplied version of lapack.  If that search failed, then try
  # to find MKL on the local system.

  if( NOT lapack_FOUND )
    if( DEFINED ENV{MKLROOT} )
      message( STATUS "Looking for lapack (MKL)...")
      # CMake uses the 'Intel10_64lp' enum to indicate MKL. For details see the
      # cmake documentation for FindBLAS.
      set( BLA_VENDOR "Intel10_64lp" )
      find_package( BLAS QUIET )

      # If we link statically, we notice that the mkl library dependencies are
      # cyclic and FindBLAS and FindLAPACK will fail.  If this is the case, but
      # we still found all the important libraries, set BLAS_FOUND=TRUE and
      # finish setting up the MKL libraries as a valid TPL for blas/lapack.
      if( NOT BLAS_FOUND AND
          BLAS_iomp5_LIBRARY AND
          BLAS_mkl_core_LIBRARY AND
          BLAS_mkl_intel_thread_LIBRARY AND
          BLAS_mkl_intel_lp64_LIBRARY )
        set( BLAS_FOUND TRUE )
      endif()

      if( "${BLAS_mkl_core_LIBRARY}" MATCHES "libmkl_core.a" )
        set( MKL_LIBRARY_TYPE "STATIC" )
      else()
        set( MKL_LIBRARY_TYPE "SHARED" )
      endif()

      if( BLAS_FOUND )
        set( LAPACK_FOUND TRUE CACHE BOOL "lapack (MKL) found?" FORCE)
        set( lapack_FOUND TRUE CACHE BOOL "lapack (MKL) found?" FORCE)
        set( lapack_DIR "$ENV{MKLROOT}" CACHE PATH "MKLROOT PATH?" FORCE)
        set( lapack_flavor "mkl")
        add_library( lapack ${MKL_LIBRARY_TYPE} IMPORTED)
        add_library( blas   ${MKL_LIBRARY_TYPE} IMPORTED)
        add_library( blas::mkl_thread  ${MKL_LIBRARY_TYPE} IMPORTED)
        add_library( blas::mkl_core    ${MKL_LIBRARY_TYPE} IMPORTED)
        set_target_properties( blas::mkl_thread PROPERTIES
          IMPORTED_LOCATION                 "${BLAS_mkl_intel_thread_LIBRARY}"
          IMPORTED_LINK_INTERFACE_LANGUAGES "C"
          IMPORTED_LINK_INTERFACE_MULTIPLICITY 20 )
        set_target_properties( blas::mkl_core PROPERTIES
          IMPORTED_LOCATION                 "${BLAS_mkl_core_LIBRARY}"
          IMPORTED_LINK_INTERFACE_LANGUAGES "C"
          IMPORTED_LINK_INTERFACE_LIBRARIES blas::mkl_thread
          IMPORTED_LINK_INTERFACE_MULTIPLICITY 20 )
        set_target_properties( blas PROPERTIES
          IMPORTED_LOCATION                 "${BLAS_mkl_intel_lp64_LIBRARY}"
          IMPORTED_LINK_INTERFACE_LANGUAGES "C"
          IMPORTED_LINK_INTERFACE_LIBRARIES "-Wl,--start-group;${BLAS_mkl_core_LIBRARY};${BLAS_mkl_intel_thread_LIBRARY};-Wl,--end-group"
          IMPORTED_LINK_INTERFACE_MULTIPLICITY 20)
        set_target_properties( lapack PROPERTIES
          IMPORTED_LOCATION                 "${BLAS_mkl_intel_lp64_LIBRARY}"
          IMPORTED_LINK_INTERFACE_LANGUAGES "C"
          IMPORTED_LINK_INTERFACE_LIBRARIES blas
          IMPORTED_LINK_INTERFACE_MULTIPLICITY 20)
        message(STATUS "Looking for lapack (MKL)...found ${BLAS_mkl_intel_lp64_LIBRARY}")
      else()
        message(STATUS "Looking for lapack (MKL)...NOTFOUND")
      endif()

    endif()
  endif()

  # If the above searches for LAPACK failed, then try to find OpenBlas on the
  # local system.

  if( NOT lapack_FOUND )
      message( STATUS "Looking for lapack (OpenBLAS)...")
      # CMake uses the 'OpenBLAS' enum to help the FindBLAS.cmake macro. For
      # details see the cmake documentation for FindBLAS.
      set( BLA_VENDOR "OpenBLAS" )
      find_package( BLAS QUIET )

      if( BLAS_FOUND )
        set( LAPACK_FOUND TRUE CACHE BOOL "lapack (OpenBlas) found?")
        set( lapack_FOUND TRUE CACHE BOOL "lapack (OpenBlas) found?")
        set( lapack_flavor "openblas")
        add_library( lapack SHARED IMPORTED)
        add_library( blas   SHARED IMPORTED)
        set_target_properties( blas PROPERTIES
          IMPORTED_LOCATION                 "${BLAS_openblas_LIBRARY}"
          IMPORTED_LINK_INTERFACE_LANGUAGES "C" )
        set_target_properties( lapack PROPERTIES
          IMPORTED_LOCATION                 "${BLAS_openblas_LIBRARY}"
          IMPORTED_LINK_INTERFACE_LANGUAGES "C" )
        message(STATUS "Looking for lapack (OpenBLAS)...found ${BLAS_openblas_LIBRARY}")
      else()
        message(STATUS "Looking for lapack (OpenBLAS)...NOTFOUND")
      endif()

  endif()

  # If the above searches for LAPACK failed, then try to find netlib-lapack and
  # netlib-blas on the local system (without the cmake config files).

  if( NOT lapack_FOUND )
      message( STATUS "Looking for lapack (no cmake config files)...")
      find_package( BLAS QUIET )

      if( BLAS_FOUND )
        find_package( LAPACK QUIET)
        set( lapack_FOUND TRUE )
        add_library( lapack SHARED IMPORTED)
        add_library( blas   SHARED IMPORTED)
        set_target_properties( blas PROPERTIES
          IMPORTED_LOCATION                 "${BLAS_blas_LIBRARY}"
          IMPORTED_LINK_INTERFACE_LANGUAGES "C" )
        set_target_properties( lapack PROPERTIES
          IMPORTED_LOCATION                 "${LAPACK_lapack_LIBRARY}"
          IMPORTED_LINK_INTERFACE_LANGUAGES "C" )
        message(STATUS "Looking for lapack(no cmake config)...found ${LAPACK_lapack_LIBRARY}")
      else()
        message(STATUS "Looking for lapack(no cmake config)...NOTFOUND")
      endif()

  endif()

  # Export BLAS/LAPACK target information to draco-config.cmake

  if( lapack_FOUND )
      # Choose items for props list via:
      # include(print_target_properties)
      # print_targets_properties("blas;lapack")
      set( props
        GNUtoMS
        IMPORTED
        IMPORTED_CONFIGURATIONS
        IMPORTED_IMPLIB
        IMPORTED_IMPLIB_DEBUG
        IMPORTED_LINK_INTERFACE_LANGUAGES
        IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG
        IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE
        IMPORTED_LINK_INTERFACE_LIBRARIES_RELWITHDEBINFO
        IMPORTED_LOCATION
        IMPORTED_LOCATION_DEBUG
        IMPORTED_LOCATION_RELEASE
        IMPORTED_LOCATION_RELWITHDEBINFO
        IMPORTED_SONAME_DEBUG
        IMPORTED_SONAME_RELEASE
        IMPORTED_SONAME_RELWITHDEBINFO
        INTERFACE_INCLUDE_DIRECTORIES
        INTERFACE_LINK_LIBRARIES
        POSITION_INDEPENDENT_CODE
        )
      if( "${lapack_flavor}" STREQUAL "mkl" )
        save_vendor_imported_library_to_draco_config(
          "blas::mkl_thread;blas::mkl_core" "${props}" )
      endif()
      save_vendor_imported_library_to_draco_config( "lapack;blas" "${props}" )
  endif()

endmacro()

#------------------------------------------------------------------------------
# Helper macros for CUDA/Unix
#
# Processes FindCUDA.cmake from the standard CMake Module location.
# This standard file establishes many macros and variables used by the
# build system for compiling CUDA code.
# (try 'cmake --help-module FindCUDA' for details)
#
# Override the FindCUDA defaults in a way that is standardized for
# Draco and Draco clients.
#
# Provided macros:
#    cuda_add_library(    target )
#    cuda_add_executable( target )
#    ...
#
# Provided variables:
#    CUDA_FOUND
#    CUDA_PROPAGATE_HOST_FLAGS
#    CUDA_NVCC_FLAGS
#    CUDA_SDK_ROOT_DIR
#    CUDA_VERBOSE_BUILD
#    CUDA_TOOLKIT_ROOT_DIR
#    CUDA_BUILD_CUBIN
#    CUDA_BUILD_EMULATION
#    ...
#------------------------------------------------------------------------------
macro( setupCudaEnv )

  if( NOT DEFINED USE_CUDA )
    option( USE_CUDA "If CUDA is available, should we use it?" OFF )
  endif()
  if( USE_CUDA )

    message( STATUS "Looking for CUDA..." )
    find_package( CUDA QUIET )
    set_package_properties( CUDA PROPERTIES
      DESCRIPTION "Toolkit providing tools and libraries needed for GPU applications."
      TYPE OPTIONAL
      PURPOSE "Required for building a GPU enabled application." )
    if( NOT EXISTS ${CUDA_NVCC_EXECUTABLE} )
      set( CUDA_FOUND 0 )
    endif()
    if( CUDA_FOUND )
      set( HAVE_CUDA 1 )
      option( USE_CUDA "If CUDA is available, should we use it?" ON )
      if( USE_CUDA )
        set( CUDA_PROPAGATE_HOST_FLAGS OFF CACHE BOOL "blah" FORCE)
        set( CUDA_NVCC_FLAGS "-arch=sm_21" )
        string( TOUPPER ${CMAKE_BUILD_TYPE} UC_CMAKE_BUILD_TYPE )
        if( ${UC_CMAKE_BUILD_TYPE} MATCHES DEBUG )
          set( CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -g -G" )
        endif()
        set( cudalibs ${CUDA_CUDART_LIBRARY} )
        set( DRACO_LIBRARY_TYPE "STATIC" CACHE STRING
          "static or shared (dll) libraries" FORCE )
      endif()
      message( STATUS "Looking for CUDA......found ${CUDA_NVCC_EXECUTABLE}" )
    else()
      message( STATUS "Looking for CUDA......not found" )
      set( USE_CUDA OFF CACHE BOOL "cuda" FORCE )
    endif()
    mark_as_advanced(
      CUDA_SDK_ROOT_DIR
      CUDA_VERBOSE_BUILD
      CUDA_TOOLKIT_ROOT_DIR
      CUDA_BUILD_CUBIN
      CUDA_BUILD_EMULATION
      CUDA_HOST_COMPILER
      )
  endif()

endmacro()

#------------------------------------------------------------------------------
# Setup QT (any)
#------------------------------------------------------------------------------
macro( setupQt )
  message( STATUS "Looking for Qt SDK..." )

  # The CMake package information should be found in
  # $QTDIR/lib/cmake/Qt5Widgets/Qt5WidgetsConfig.cmake.  On CCS Linux
  # machines, QTDIR is set when loading the qt module
  # (QTDIR=/ccs/codes/radtran/vendors/Qt53/5.3/gcc_64):
  if( "${QTDIR}notset" STREQUAL "notset" AND EXISTS "$ENV{QTDIR}" )
    set( QTDIR $ENV{QTDIR} CACHE PATH "This path should include /lib/cmake/Qt5Widgets" )
  endif()
  set( CMAKE_PREFIX_PATH_QT "$ENV{QTDIR}/lib/cmake/Qt5Widgets" )

  if( NOT EXISTS ${CMAKE_PREFIX_PATH_QT}/Qt5WidgetsConfig.cmake )
    # message( FATAL_ERROR "Could not find cQt cmake macros.  Try
    # setting CMAKE_PREFIX_PATH_QT to the path that contains
    # Qt5WidgetsConfig.cmake" )
    message( STATUS "Looking for Qt SDK...not found." )
  else()
    file( TO_CMAKE_PATH "${CMAKE_PREFIX_PATH_QT}" CMAKE_PREFIX_PATH_QT )
    list( APPEND CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH_QT}" )
    find_package(Qt5Widgets)
    find_package(Qt5Core)
    get_target_property(QtCore_location Qt5::Core LOCATION)
    if( Qt5Widgets_FOUND )
      set( QT_FOUND 1 )
      # Instruct CMake to run moc automatically when needed (only for
      # subdirectories that need Qt)
      # set(CMAKE_AUTOMOC ON)
      message( STATUS "Looking for Qt SDK...${QTDIR}." )
    else()
      set( QT_FOUND "QT-NOTFOUND" )
      message( STATUS "Looking for Qt SDK...not found." )
    endif()
  endif()

  if( QT_FOUND )
    mark_as_advanced( Qt5Core_DIR Qt5Gui_DIR Qt5Gui_EGL_LIBRARY
      Qt5Widgets_DIR )
  endif()
endmacro()

#------------------------------------------------------------------------------
# save_vendor_imported_library_to_draco_config
#
# Save imported target information to the installed file
# "draco-config.cmake"
#
# Append the variable Draco_EXPORT_TARGET_PROPERTIES with enough
# information to use the imported library from Jayenne or
# Capsaicin.  The string  Draco_EXPORT_TARGET_PROPERTIES is written
# to lib/cmake/draco/draco-config.cmake during the 'install' phase
# of the build.
# - Examine each imported target specified in the the list 'targets'
#   and generate an 'import_library' for the imported target.
# - Query each target for each of the properites listed in
#   'target_properties'
# - If the target has the requested property, save this value to
#   Draco_EXPORT_TARGET_PROPERTIES.
#------------------------------------------------------------------------------
macro( save_vendor_imported_library_to_draco_config targets target_properties )
  foreach( tgt ${targets} )
    # The target's TYPE may not be set correctly, set it manually.
    get_target_property( fplib ${tgt} IMPORTED_LOCATION )
    get_target_property( fplibr ${tgt} IMPORTED_LOCATION_RELEASE )
    get_target_property( fplibd ${tgt} IMPORTED_LOCATION_DEBUG )
    set( fplib "${fplib} ${fplibr} ${fplibd}" )
    if( "${fplib}" MATCHES "[.]so"  OR
        "${fplib}" MATCHES "[.]dll" OR
        "${fplib}" MATCHES "[.]dylib" )
      set( library_type SHARED )
    else()
      set( library_type STATIC )
      #set( library_type UNKNOWN )
    endif()

    set( tmp "add_library( ${tgt} ${library_type} IMPORTED )" )
    set( tmp "${tmp}
set_target_properties( ${tgt} PROPERTIES")
    if(TARGET ${tgt})
      foreach( prop ${target_properties} )
        get_property(v TARGET ${tgt} PROPERTY ${prop})
        #get_property(d TARGET ${tgt} PROPERTY ${prop} DEFINED)
        get_property(s TARGET ${tgt} PROPERTY ${prop} SET)
        if( s )
          set( tmp "${tmp}
   ${prop} \"${v}\"")
        endif()
      endforeach()
      set( tmp "${tmp} )
")
    else()
      message(FATAL_ERROR "There is no target named '${tgt}'")
    endif()

    # Append to Draco_EXPORT_TARGET_PROPERTIES the commands needed to
    # import the library and set all of its properties.
    set( Draco_EXPORT_TARGET_PROPERTIES
      "${Draco_EXPORT_TARGET_PROPERTIES}
${tmp}")
    # message("${tmp}")
  endforeach()
  unset( fplib )
  unset( library_type )
endmacro()

#------------------------------------------------------------------------------
# Setup GSL (any)
#------------------------------------------------------------------------------
macro( setupGSL )
  if( NOT TARGET GSL::gsl )

    message( STATUS "Looking for GSL..." )

    # If gsl-config is in the PATH, query the value for GSL_ROOT_DIR
    # This bit of logic is needed on Cielo/Cielito because gsl is not in
    # a system location (it is provided by a module)
    if( "$ENV{GSL_ROOT_DIR}x" STREQUAL "x" AND "${GSL_ROOT_DIR}x" STREQUAL "x")
      find_program( GSL_CONFIG gsl-config )
      if( EXISTS "${GSL_CONFIG}" )
        exec_program( "${GSL_CONFIG}"
          ARGS --prefix
          OUTPUT_VARIABLE GSL_ROOT_DIR )
      endif()
    endif()

    find_package( GSL QUIET REQUIRED )
    if( GSL_FOUND )
      message( STATUS "Looking for GSL.......found ${GSL_LIBRARY}" )
      mark_as_advanced( GSL_CONFIG_EXECUTABLE )

      # Export GSL target information to draco-config.cmake
      # Choose items for props list via:
      # include(print_target_properties)
      # print_targets_properties("GSL::gsl;GSL::gslcblas")
      set( props
        GNUtoMS
        IMPORTED_CONFIGURATIONS
        IMPORTED_IMPLIB
        IMPORTED_IMPLIB_DEBUG
        IMPORTED_LINK_INTERFACE_LANGUAGES
        INTERFACE_LINK_LIBRARIES
        IMPORTED_LOCATION_DEBUG
        IMPORTED_LOCATION_RELEASE
        IMPORTED
        INTERFACE_INCLUDE_DIRECTORIES
        POSITION_INDEPENDENT_CODE
        IMPORTED_LOCATION )
      save_vendor_imported_library_to_draco_config(
        "GSL::gsl;GSL::gslcblas" "${props}" )
    else()
      message( STATUS "Looking for GSL.......not found" )
    endif()

    #=============================================================================
    # Include some information that can be printed by the build system.
    set_package_properties( GSL PROPERTIES
      DESCRIPTION "Gnu Scientific Library"
      URL "www.gnu.org/software/gsl"
      PURPOSE "The GNU Scientific Library (GSL) is a numerical library for C and C++ programmers."
      )

  endif()

endmacro()

#------------------------------------------------------------------------------
# Setup ParMETIS (any)
#------------------------------------------------------------------------------
macro( setupParMETIS )

  if( NOT TARGET METIS::metis )
    message( STATUS "Looking for METIS..." )

    find_package( METIS QUIET )
    if( METIS_FOUND )
      message( STATUS "Looking for METIS.....found ${METIS_LIBRARY}" )

      # Export METIS target information to draco-config.cmake
      # Choose items for props list via:
      # include(print_target_properties)
      # print_targets_properties("METIS::metis")
      set( props
        GNUtoMS
        IMPORTED_CONFIGURATIONS
        IMPORTED_IMPLIB
        IMPORTED_IMPLIB_DEBUG
        IMPORTED_LINK_INTERFACE_LANGUAGES
        INTERFACE_LINK_LIBRARIES
        IMPORTED_LOCATION_DEBUG
        IMPORTED_LOCATION_RELEASE
        IMPORTED
        INTERFACE_INCLUDE_DIRECTORIES
        POSITION_INDEPENDENT_CODE
        IMPORTED_LOCATION )

      save_vendor_imported_library_to_draco_config( "METIS::metis" "${props}" )
    else()
      message( STATUS "Looking for METIS.....not found" )
    endif()

    #=============================================================================
    # Include some information that can be printed by the build system.
    set_package_properties( METIS PROPERTIES
      DESCRIPTION "METIS"
      URL "http://glaros.dtc.umn.edu/gkhome/metis/metis/overview"
      PURPOSE "METIS is a set of serial programs for partitioning graphs, partitioning finite
   element meshes, and producing fill reducing orderings for sparse matrices."
      )

  endif()

  if( NOT TARGET ParMETIS::parmetis )

    message( STATUS "Looking for ParMETIS..." )

    find_package( ParMETIS QUIET )
    if( ParMETIS_FOUND )
      message( STATUS "Looking for ParMETIS..found ${ParMETIS_LIBRARY}" )

      # Export ParMETIS target information to draco-config.cmake
      # Choose items for props list via:
      # include(print_target_properties)
      # print_targets_properties("ParMETIS::parmetis")
      set( props
        GNUtoMS
        IMPORTED_CONFIGURATIONS
        IMPORTED_IMPLIB
        IMPORTED_IMPLIB_DEBUG
        IMPORTED_LINK_INTERFACE_LANGUAGES
        INTERFACE_LINK_LIBRARIES
        IMPORTED_LOCATION_DEBUG
        IMPORTED_LOCATION_RELEASE
        IMPORTED
        INTERFACE_INCLUDE_DIRECTORIES
        POSITION_INDEPENDENT_CODE
        IMPORTED_LOCATION )

      save_vendor_imported_library_to_draco_config(
        "ParMETIS::parmetis" "${props}" )
    else()
      message( STATUS "Looking for ParMETIS..not found" )
    endif()

    #=============================================================================
    # Include some information that can be printed by the build system.
    set_package_properties( ParMETIS PROPERTIES
      DESCRIPTION "MPI Parallel METIS"
      URL "http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview"
      PURPOSE "ParMETIS is an MPI-based parallel library that implements a variety of algorithms for partitioning
   unstructured graphs, meshes, and for computing fill-reducing orderings of sparse matrices."
      )

  endif()

endmacro()

#------------------------------------------------------------------------------
# Setup SuperLU_DIST (any)
#------------------------------------------------------------------------------
macro( setupSuperLU_DIST )

  if( NOT TARGET SuperLU_DIST::superludist )
    message( STATUS "Looking for SuperLU_DIST..." )

    find_package( SuperLU_DIST QUIET )
    if( SuperLU_DIST_FOUND )
      message( STATUS "Looking for SuperLU_DIST.....found ${SuperLU_DIST_LIBRARY}" )

      # Export SuperLU_DIST target information to draco-config.cmake
      # Choose items for props list via:
      # include(print_target_properties)
      # print_targets_properties("SuperLU_DIST::superludist")
      set( props
        GNUtoMS
        IMPORTED_CONFIGURATIONS
        IMPORTED_IMPLIB
        IMPORTED_IMPLIB_DEBUG
        IMPORTED_LINK_INTERFACE_LANGUAGES
        INTERFACE_LINK_LIBRARIES
        IMPORTED_LOCATION_DEBUG
        IMPORTED_LOCATION_RELEASE
        IMPORTED
        INTERFACE_INCLUDE_DIRECTORIES
        POSITION_INDEPENDENT_CODE
        IMPORTED_LOCATION )

      save_vendor_imported_library_to_draco_config( "SuperLU_DIST::superludist" "${props}" )
    else()
      message( STATUS "Looking for SuperLU_DIST.....not found" )
    endif()

    #=============================================================================
    # Include some information that can be printed by the build system.
    set_package_properties( SuperLU_DIST PROPERTIES
      DESCRIPTION "SuperLU_DIST"
      URL " http://crd-legacy.lbl.gov/~xiaoye/SuperLU/"
      PURPOSE "SuperLU is a general purpose library for the direct solution of large, sparse,
   nonsymmetric systems of linear equations on high performance machines."
      )

  endif()

endmacro()

#------------------------------------------------------------------------------
# Helper macros for setup_global_libraries()
#------------------------------------------------------------------------------
macro( SetupVendorLibrariesUnix )

  # GSL, METIS, ParMETIS and SuperLU_DIST ------------------------------------
  setupGSL()
  setupParMETIS()
  setupSuperLU_DIST()

  # Random123 ----------------------------------------------------------------
  message( STATUS "Looking for Random123...")
  find_package( Random123 REQUIRED QUIET )
  if( RANDOM123_FOUND )
    message( STATUS "Looking for Random123.found ${RANDOM123_INCLUDE_DIR}")
  else()
    message( STATUS "Looking for Random123.not found")
  endif()

  # Grace ------------------------------------------------------------------
  message( STATUS "Looking for Grace...")
  find_package( Grace QUIET )
  set_package_properties( Grace PROPERTIES
    DESCRIPTION "A WYSIWYG 2D plotting tool."
    TYPE OPTIONAL
    PURPOSE "Required for building the plot2D component."
    )
  if( Grace_FOUND )
    message( STATUS "Looking for Grace.....found ${Grace_EXECUTABLE}")
  else()
    message( STATUS "Looking for Grace.....not found")
  endif()

  # CUDA ------------------------------------------------------------------
  setupCudaEnv()

  # PYTHON ----------------------------------------------------------------

  message( STATUS "Looking for Python...." )
  find_package(PythonInterp QUIET)
  #  PYTHONINTERP_FOUND - Was the Python executable found
  #  PYTHON_EXECUTABLE  - path to the Python interpreter
  set_package_properties( PythonInterp PROPERTIES
    DESCRIPTION "Python interpreter"
    TYPE OPTIONAL
    PURPOSE "Required for running the fpe_trap tests."
    )
  if( PYTHONINTERP_FOUND )
    # if( ${PYTHON_VERSION_STRING} GREATER 3.9.9 )
    #   message( FATAL_ERROR
    #     "Looking for Python....found version"
    #     "${PYTHON_VERSION_STRING}, but Draco requires version < 4.0." )
    # endif()
    message( STATUS "Looking for Python....found ${PYTHON_EXECUTABLE}" )
  else()
    message( STATUS "Looking for Python....not found" )
  endif()

  # Qt -----------------------------------------------------------------------
  setupQt()

endmacro()

##---------------------------------------------------------------------------##
## Vendors for building on Windows-based platforms.
##---------------------------------------------------------------------------##

macro( SetupVendorLibrariesWindows )

  # GSL ---------------------------------------------------------------------
  setupGSL()
  setupParMETIS()

  # Random123 ---------------------------------------------------------------
  message( STATUS "Looking for Random123...")
  find_package( Random123 REQUIRED QUIET )
  if( RANDOM123_FOUND )
    message( STATUS "Looking for Random123.found ${RANDOM123_INCLUDE_DIR}")
  else()
    message( STATUS "Looking for Random123.not found")
  endif()

  # PYTHON ----------------------------------------------------------------
  find_package(PythonInterp QUIET)
  #  PYTHONINTERP_FOUND - Was the Python executable found
  #  PYTHON_EXECUTABLE  - path to the Python interpreter
  set_package_properties( PythonInterp PROPERTIES
    DESCRIPTION "Python interpreter"
    TYPE OPTIONAL
    PURPOSE "Required for running the fpe_trap tests."
    )

  # Qt -----------------------------------------------------------------------
  setupQt()

endmacro()

#------------------------------------------------------------------------------
# Helper macros for setup_global_libraries()
# Assign here the library version to be used.
#------------------------------------------------------------------------------
macro( setVendorVersionDefaults )
  #Set the preferred search directories(ROOT)

  #Check that VENDOR_DIR is defined as a cache variable or as an
  #environment variable. If defined as both then take the
  #environment variable.

  # See if VENDOR_DIR is set.  Try some defaults if it is not set.
  if( NOT EXISTS "${VENDOR_DIR}" AND IS_DIRECTORY "$ENV{VENDOR_DIR}" )
    set( VENDOR_DIR $ENV{VENDOR_DIR} )
  endif()
  # If needed, try some obvious places.
  if( NOT EXISTS "${VENDOR_DIR}" )
    if( IS_DIRECTORY "/ccs/codes/radtran/vendors/Linux64" )
      set( VENDOR_DIR "/ccs/codes/radtran/vendors/Linux64" )
    endif()
    if( IS_DIRECTORY /usr/projects/draco/vendors )
      set( VENDOR_DIR /usr/projects/draco/vendors )
    endif()
    if( IS_DIRECTORY c:/vendors )
      set( VENDOR_DIR c:/vendors )
    endif()
  endif()
  # Cache the result
  if( IS_DIRECTORY "${VENDOR_DIR}")
    set( VENDOR_DIR ${VENDOR_DIR} CACHE PATH
      "Root directory where CCS-2 3rd party libraries are located."
      FORCE )
  else()
    message( "
WARNING: VENDOR_DIR not defined locally or in user environment,
individual vendor directories should be defined." )
  endif()

  # Import environment variables related to vendors
  # 1. Use command line variables (-DLAPACK_LIB_DIR=<path>
  # 2. Use environment variables ($ENV{LAPACK_LIB_DIR}=<path>)
  # 3. Try to find vendor in $VENDOR_DIR
  # 4. Don't set anything and let the user set a value in the cache
  #    after failed 1st configure attempt.
  if( NOT LAPACK_LIB_DIR AND IS_DIRECTORY $ENV{LAPACK_LIB_DIR} )
    set( LAPACK_LIB_DIR $ENV{LAPACK_LIB_DIR} )
    set( LAPACK_INC_DIR $ENV{LAPACK_INC_DIR} )
  endif()
  if( NOT LAPACK_LIB_DIR AND IS_DIRECTORY ${VENDOR_DIR}/lapack-3.4.2/lib )
    set( LAPACK_LIB_DIR "${VENDOR_DIR}/lapack-3.4.2/lib" )
    set( LAPACK_INC_DIR "${VENDOR_DIR}/lapack-3.4.2/include" )
  endif()
  # if( NOT LAPACK_LIB_DIR AND IS_DIRECTORY ${VENDOR_DIR}/clapack/lib )
  # set( LAPACK_LIB_DIR "${VENDOR_DIR}/clapack/lib" )
  # set( LAPACK_INC_DIR "${VENDOR_DIR}/clapack/include" )
  # endif()

  if( NOT GSL_LIB_DIR )
    if( IS_DIRECTORY $ENV{GSL_LIB_DIR}  )
      set( GSL_LIB_DIR $ENV{GSL_LIB_DIR} )
      set( GSL_INC_DIR $ENV{GSL_INC_DIR} )
    elseif( IS_DIRECTORY ${VENDOR_DIR}/gsl/lib )
      set( GSL_LIB_DIR "${VENDOR_DIR}/gsl/lib" )
      set( GSL_INC_DIR "${VENDOR_DIR}/gsl/include" )
    endif()
  endif()

  if( NOT ParMETIS_ROOT_DIR )
    if( IS_DIRECTORY $ENV{ParMETIS_ROOT_DIR}  )
      set( ParMETIS_ROOT_DIR $ENV{ParMETIS_ROOT_DIR} )
    endif()
  endif()

  if( NOT RANDOM123_INC_DIR AND IS_DIRECTORY $ENV{RANDOM123_INC_DIR}  )
    set( RANDOM123_INC_DIR $ENV{RANDOM123_INC_DIR} )
  endif()
  if( NOT RANDOM123_INC_DIR AND IS_DIRECTORY ${VENDOR_DIR}/Random123-1.08/include )
    set( RANDOM123_INC_DIR "${VENDOR_DIR}/Random123-1.08/include" )
  endif()

  set_package_properties( BLAS PROPERTIES
    DESCRIPTION "Basic Linear Algebra Subprograms"
    TYPE OPTIONAL
    PURPOSE "Required for building the lapack_wrap component."
    )
  set_package_properties( lapack PROPERTIES
    DESCRIPTION "Linear Algebra PACKage"
    TYPE OPTIONAL
    PURPOSE "Required for building the lapack_wrap component."
    )
  set_package_properties( Random123 PROPERTIES
    URL "http://www.deshawresearch.com/resources_random123.html"
    DESCRIPTION "a library of counter-based random number generators
"
    TYPE REQUIRED
    PURPOSE "Required for building rng component."
    )

endmacro()

#------------------------------------------------------------------------------
# This macro should contain all the system libraries which are
# required to link the main objects.
#------------------------------------------------------------------------------
macro( setupVendorLibraries )

  message( "\nVendor Setup:\n")

  #
  # General settings
  #
  setVendorVersionDefaults()
  if( NOT TARGET blas )
    setupLAPACKLibrariesUnix()
  endif()

  # System specific settings
  if ( UNIX )

    if( NOT MPI_CXX_COMPILER )
      setupMPILibrariesUnix()
    endif()
    setupVendorLibrariesUnix()

  elseif( WIN32 )

    if( NOT MPI_CXX_COMPILER )
      setupMPILibrariesWindows()
    endif()
    setupVendorLibrariesWindows()

  else()
    message( FATAL_ERROR "
I don't know how to setup global (vendor) libraries for this platform.
WIN32=0; UNIX=0; CMAKE_SYSTEM=${CMAKE_SYSTEM};
CMAKE_SYSTEM_NAME=${CMAKE_SYSTEM_NAME}" )
  endif()

endmacro()

#----------------------------------------------------------------------#
# End vendor_libraries.cmake
#----------------------------------------------------------------------#
