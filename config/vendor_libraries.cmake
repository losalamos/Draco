#-----------------------------*-cmake-*----------------------------------------#
# file   config/vendor_libraries.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 June 6
# brief  Setup Vendors
# note   Copyright (C) 2010-2013 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$ 
#------------------------------------------------------------------------------#

#
# Look for any libraries which are required at the toplevel.
# 

include( FeatureSummary )
include( setupMPI ) # defines the macros setupMPILibrariesUnix|Windows

##---------------------------------------------------------------------------##
## Items to consider...
##---------------------------------------------------------------------------##

#
# Getting .so requirements?
#
#include( GetPrerequisites )
#get_prerequisites( library LIST_OF_REQS 1 0 <exepath> <dirs>)
#list_prerequisites(<target> [<recurse> [<exclude_system> [<verbose>]]])

#
# Install all required system libraries?
#
# include( InstallRequiredSystemLibraries)

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

  message( STATUS "Looking for lapack...")
  set( lapack_FOUND FALSE )
  # Use LAPACK_LIB_DIR, if the user set it, to help find LAPACK. 
  foreach( version 3.4.1 3.4.2 3.5.0 )   
    if( EXISTS  ${LAPACK_LIB_DIR}/cmake/lapack-${version} )
      list( APPEND CMAKE_PREFIX_PATH ${LAPACK_LIB_DIR}/cmake/lapack-${version} )
      find_package( lapack CONFIG )
    endif()
  endforeach()
  if( lapack_FOUND )
    foreach( config NOCONFIG DEBUG RELEASE )
      get_target_property(tmp lapack IMPORTED_LOCATION_${config} )
      if( EXISTS ${tmp} )
        set( lapack_FOUND TRUE )
      endif()
    endforeach()
    message( STATUS "Looking for lapack....found ${LAPACK_LIB_DIR}")
    set( lapack_FOUND ${lapack_FOUND} CACHE BOOL "Did we find LAPACK." FORCE )
  else()
    message( STATUS "Looking for lapack....not found")
  endif()

  mark_as_advanced( lapack_DIR lapack_FOUND )

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

  message( STATUS "Looking for CUDA..." )
  # special code for CT/CI
  if( "${CMAKE_SYSTEM_PROCESSOR}notset" STREQUAL "notset" AND 
      ${CMAKE_SYSTEM_NAME} MATCHES "Catamount")
    set( CMAKE_SYSTEM_PROCESSOR "x86_64" CACHE STRING 
      "For unix, this value is set from uname -p." FORCE)
  endif()
  find_package( CUDA QUIET )
  set_package_properties( CUDA PROPERTIES
    DESCRIPTION "Toolkit providing tools and libraries needed for GPU applications."
    TYPE OPTIONAL
    PURPOSE "Required for bulding a GPU enabled application." )
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
  endif()
  mark_as_advanced( 
    CUDA_SDK_ROOT_DIR 
    CUDA_VERBOSE_BUILD
    CUDA_TOOLKIT_ROOT_DIR 
    CUDA_BUILD_CUBIN
    CUDA_BUILD_EMULATION
    )

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
endmacro()

#------------------------------------------------------------------------------
# Setup GSL (any)
#------------------------------------------------------------------------------
macro( setupGSL )
  message( STATUS "Looking for GSL..." )
  find_package( GSL QUIET REQUIRED )
  if( GSL_FOUND )
    message( STATUS "Looking for GSL.......found ${GSL_LIBRARY}" )

    # Create an entry in draco-config.cmake for the gsl libs
    get_target_property( gslimploc      gsl::gsl      IMPORTED_LOCATION )
    get_target_property( gslcblasimploc gsl::gslcblas IMPORTED_LOCATION )
    
    # If this platform doesn't support shared libraries (e.g. cross
    # compiling), assume static. This suppresses cmake (3.0.0) warnings
    # of the form: 
    #     "ADD_LIBRARY called with MODULE option but the target platform
    #     does not support dynamic linking.  Building a STATIC library
    #     instead.  This may lead to problems."
  if( TARGET_SUPPORTS_SHARED_LIBS )
    set( library_type UNKNOWN )
  else()
    set( library_type STATIC )
  endif()

    set( Draco_EXPORT_TARGET_PROPERTIES 
      "${Draco_EXPORT_TARGET_PROPERTIES}
add_library( gsl::gsl ${library_type} IMPORTED )
set_target_properties( gsl::gsl PROPERTIES
    IMPORTED_LINK_INTERFACE_LANGUAGES \"C\"
    IMPORTED_LINK_INTERFACE_LIBRARIES \"gsl::gslcblas\"
    IMPORTED_LOCATION                 \"${gslimploc}\" 
)

add_library( gsl::gslcblas ${library_type} IMPORTED )
set_target_properties( gsl::gslcblas PROPERTIES
    IMPORTED_LINK_INTERFACE_LANGUAGES \"C\"
    IMPORTED_LOCATION                 \"${gslcblasimploc}\" 
)
")    

  unset(library_type)

else()
  message( STATUS "Looking for GSL.......not found" )
endif()
endmacro()

#------------------------------------------------------------------------------
# Helper macros for setup_global_libraries()
#------------------------------------------------------------------------------
macro( SetupVendorLibrariesUnix )

  # GSL ----------------------------------------------------------------------
  setupGSL()
  
  # Random123 ----------------------------------------------------------------
  message( STATUS "Looking for Random123...")
  find_package( Random123 REQUIRED QUIET )
  if( RANDOM123_FOUND )
    message( STATUS "Looking for Random123.found ${RANDOM123_INCLUDE_DIR}")
  else()
    message( STATUS "Looking for Random123.not found")
  endif()

  # GRACE ------------------------------------------------------------------
  find_package( Grace QUIET )
  set_package_properties( Grace PROPERTIES
    DESCRIPTION "A WYSIWYG 2D plotting tool."
    TYPE OPTIONAL
    PURPOSE "Required for bulding the plot2D component."
    )

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
    set( VENDOR_DIR $ENV{VENDOR_DIR} CACHE PATH
      "Root directory where CCS-2 3rd party libraries are located." )
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
    elseif( IS_DIRECTORY c:/a/vendors )
      set( VENDOR_DIR c:/a/vendors )
    endif()
  endif()
  # Cache the result
  if( IS_DIRECTORY "${VENDOR_DIR}")
    set( VENDOR_DIR $ENV{VENDOR_DIR} CACHE PATH
      "Root directory where CCS-2 3rd party libraries are located." )
  else( IS_DIRECTORY "${VENDOR_DIR}")
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

  if( NOT RANDOM123_INC_DIR AND IS_DIRECTORY $ENV{RANDOM123_INC_DIR}  )
    set( RANDOM123_INC_DIR $ENV{RANDOM123_INC_DIR} )
  endif()
  if( NOT RANDOM123_INC_DIR AND IS_DIRECTORY ${VENDOR_DIR}/Random123-1.08/include )
    set( RANDOM123_INC_DIR "${VENDOR_DIR}/Random123-1.08/include" )
  endif()
  
  set_package_properties( MPI PROPERTIES
    URL "http://www.open-mpi.org/"
    DESCRIPTION "A High Performance Message Passing Library"
    TYPE RECOMMENDED
    PURPOSE "If not available, all Draco components will be built as scalar applications."
    )
  set_package_properties( BLAS PROPERTIES
    DESCRIPTION "Basic Linear Algebra Subprograms"
    TYPE OPTIONAL
    PURPOSE "Required for bulding the lapack_wrap component." 
    )
  set_package_properties( lapack PROPERTIES
    DESCRIPTION "Linear Algebra PACKage"
    TYPE OPTIONAL
    PURPOSE "Required for bulding the lapack_wrap component." 
    )     
  # set_package_properties( GSL PROPERTIES
  #    URL "http://www.gnu.org/s/gsl/"
  #    DESCRIPTION "GNU Scientific Library"
  #    TYPE REQUIRED
  #    PURPOSE "Required for bulding quadrature and rng components."
  #    )  
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

  # System specific settings
  if ( UNIX )
    
    if( NOT MPI_SETUP_DONE )
      setupMPILibrariesUnix()
    endif()
    setupLAPACKLibrariesUnix()
    setupVendorLibrariesUnix()
    
  elseif( WIN32 )
    
    setupMPILibrariesWindows()
    setupLAPACKLibrariesUnix()
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
