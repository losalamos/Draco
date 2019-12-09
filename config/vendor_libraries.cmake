#-----------------------------*-cmake-*----------------------------------------#
# file   config/vendor_libraries.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 June 6
# brief  Look for any libraries which are required at the top level.
# note   Copyright (C) 2016-2019 Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

include_guard(GLOBAL)
include( FeatureSummary )
include( setupMPI ) # defines the macros setupMPILibrariesUnix|Windows

#------------------------------------------------------------------------------#
# Helper macros for Python
#------------------------------------------------------------------------------#
macro( setupPython )

  message( STATUS "Looking for Python...." )
  find_package(Python QUIET REQUIRED COMPONENTS Interpreter)
  #  Python_Interpreter_FOUND - Was the Python executable found
  #  Python_EXECUTABLE  - path to the Python interpreter
  set_package_properties( PythonInterp PROPERTIES
    URL "https://www.python.org"
    DESCRIPTION "Python interpreter"
    TYPE REQUIRED
    PURPOSE "Required for running tests and accessing features that rely on matplotlib."
    )
  if( Python_Interpreter_FOUND )
    message( STATUS "Looking for Python....found ${Python_EXECUTABLE}" )
  else()
    message( STATUS "Looking for Python....not found" )
  endif()

endmacro()

#------------------------------------------------------------------------------#
# Helper macros for Random123
#
# Providers: Linux - use spack to install netlib-lapack
#                    https://github.com/spack/spack
#------------------------------------------------------------------------------#
macro( setupRandom123 )

 message( STATUS "Looking for Random123...")
  find_package( Random123 REQUIRED QUIET )
  mark_as_advanced( RANDOM123_FOUND )
  if( RANDOM123_FOUND )
    message( STATUS "Looking for Random123.found ${RANDOM123_INCLUDE_DIR}")
  else()
    message( STATUS "Looking for Random123.not found")
  endif()
  set_package_properties( Random123 PROPERTIES
    URL "http://www.deshawresearch.com/resources_random123.html"
    DESCRIPTION "a library of counter-based random number generators"
    TYPE REQUIRED
    PURPOSE "Required for building the rng component."  )
endmacro()

#------------------------------------------------------------------------------
# Helper macros for LAPACK/Unix
#
# This module sets the following variables:
# lapack_FOUND - set to true if a library implementing the LAPACK
#         interface is found
# lapack_VERSION - '3.4.1'
# provides targets: lapack, blas
#
# Providers: Linux - use spack to install netlib-lapack
#                    https://github.com/spack/spack
#            Windows - clone and build from sources
#                    https://github.com/KineticTheory/lapack-visualstudio-mingw-gfortran
#------------------------------------------------------------------------------
macro( setupLAPACKLibraries )

  # defaults
  set( lapack_url "http://www.netlib.org/lapack" )
  set( LAPACK_FOUND FALSE ) # for robustness, always do this search.

  # There are several flavors of LAPACK.
  # 1. look for netlib-lapack
  # 2. look for MKL (Intel)
  # 3. look for OpenBLAS.

  #----------------------------------------------------------------------------#
  # netlib-lapack (find_package config mode search)
  message( STATUS "Looking for lapack (netlib)...")

  #----------------------------------------------------------------------------#
  # Package-config mode: look for lapack-config.cmake in $CMAKE_PREFIX_PATH
  if( NOT TARGET lapack AND NOT "${LAPACK_FOUND}" )
    find_package( lapack CONFIG QUIET )
  endif()

  if( NOT TARGET lapack )
    message( STATUS "Looking for lapack (netlib)....not found")
  else()
    set( lapack_flavor "netlib")
    foreach( config NOCONFIG DEBUG RELEASE RELWITHDEBINFO )
      get_target_property(tmp lapack IMPORTED_LOCATION_${config} )
      if( EXISTS ${tmp} )
        set( lapack_loc ${tmp} )
        break()
      endif()
    endforeach()
    message( STATUS "Looking for lapack (netlib)....found ${lapack_loc}")

    # The above might define blas, or it might not. Double check:
    if( NOT TARGET blas )
      find_package( BLAS QUIET)
      if( BLAS_FOUND )
        add_library( blas STATIC IMPORTED)
        set_target_properties( blas PROPERTIES
          IMPORTED_LOCATION                 "${BLAS_LIBRARIES}"
          IMPORTED_LINK_INTERFACE_LANGUAGES "Fortran")
      else()
        message( FATAL_ERROR "Looking for lapack (netlib)....blas not found")
      endif()
    else()
      # ensure lapack --> blas?
      get_target_property( ilil lapack IMPORTED_LINK_INTERFACE_LIBRARIES )
      if( NOT "${ilil}" MATCHES "blas" )
        set_target_properties( lapack PROPERTIES
          IMPORTED_LINK_INTERFACE_LIBRARIES blas )
      endif()
    endif()
  endif()

  #----------------------------------------------------------------------------#
  # MKL

  # If the above search failed, then try to find MKL on the local system.
  if( NOT TARGET lapack AND NOT "${LAPACK_FOUND}" )
    if( DEFINED ENV{MKLROOT} )
      message( STATUS "Looking for lapack (MKL)...")
      # CMake uses the 'Intel10_64lp' enum to indicate MKL. For details see the
      # cmake documentation for FindBLAS.
      set( BLA_VENDOR "Intel10_64lp" )
      find_package( Threads QUIET )
      find_package( BLAS QUIET)
      find_package( LAPACK QUIET)

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

      # should we link against libmkl_gnu_thread.so or libmkl_intel_thread.so
      if( ${CMAKE_C_COMPILER_ID} MATCHES GNU )
        set(tlib "gnu")
        set(lplib "gf")
      else()
        set(tlib "intel")
        set(lplib "intel")
      endif()

      if( BLAS_FOUND )
        # set( BLAS_FOUND TRUE CACHE BOOL "lapack (MKL) found?" FORCE)
        set( LAPACK_FOUND TRUE CACHE BOOL "lapack (MKL) found?" FORCE)
        set( lapack_DIR "$ENV{MKLROOT}" CACHE PATH "MKLROOT PATH?" FORCE)
        set( lapack_flavor "mkl")
        set( lapack_url "https://software.intel.com/en-us/intel-mkl")
        add_library( lapack ${MKL_LIBRARY_TYPE} IMPORTED)
        add_library( blas   ${MKL_LIBRARY_TYPE} IMPORTED)
        add_library( blas::mkl_thread  ${MKL_LIBRARY_TYPE} IMPORTED)
        add_library( blas::mkl_core    ${MKL_LIBRARY_TYPE} IMPORTED)
        set_target_properties( blas::mkl_thread PROPERTIES
          IMPORTED_LOCATION                 "${BLAS_mkl_${tlib}_thread_LIBRARY}"
          IMPORTED_LINK_INTERFACE_LANGUAGES "C"
          IMPORTED_LINK_INTERFACE_MULTIPLICITY 20 )
        set_target_properties( blas::mkl_core PROPERTIES
          IMPORTED_LOCATION                 "${BLAS_mkl_core_LIBRARY}"
          IMPORTED_LINK_INTERFACE_LANGUAGES "C"
          IMPORTED_LINK_INTERFACE_LIBRARIES blas::mkl_thread
          IMPORTED_LINK_INTERFACE_MULTIPLICITY 20 )
        set_target_properties( blas PROPERTIES
          IMPORTED_LOCATION                 "${BLAS_mkl_${lplib}_lp64_LIBRARY}"
          IMPORTED_LINK_INTERFACE_LANGUAGES "C"
          IMPORTED_LINK_INTERFACE_LIBRARIES blas::mkl_core
#          IMPORTED_LINK_INTERFACE_LIBRARIES "-Wl,--start-group;${BLAS_mkl_core_LIBRARY};${BLAS_mkl_${tlib}_thread_LIBRARY};-Wl,--end-group"
          IMPORTED_LINK_INTERFACE_MULTIPLICITY 20)
        set_target_properties( lapack PROPERTIES
          IMPORTED_LOCATION                 "${BLAS_mkl_${lplib}_lp64_LIBRARY}"
          IMPORTED_LINK_INTERFACE_LANGUAGES "C"
          IMPORTED_LINK_INTERFACE_LIBRARIES blas
          IMPORTED_LINK_INTERFACE_MULTIPLICITY 20)
        message(STATUS "Looking for lapack (MKL)...found ${BLAS_mkl_${lplib}_lp64_LIBRARY}")
      else()
        message(STATUS "Looking for lapack (MKL)...NOTFOUND")
      endif()

    endif()
  endif()

  #----------------------------------------------------------------------------#
  # OpenBLAS

  # If the above searches for LAPACK failed, then try to find OpenBlas on the
  # local system.

  if( NOT TARGET lapack AND NOT "${LAPACK_FOUND}" )
      message( STATUS "Looking for lapack (OpenBLAS)...")
      # CMake uses the 'OpenBLAS' enum to help the FindBLAS.cmake macro. For
      # details see the cmake documentation for FindBLAS.
      set( BLA_VENDOR "OpenBLAS" )
      find_package( BLAS QUIET )

      if( BLAS_FOUND )
        set( LAPACK_FOUND TRUE CACHE BOOL "lapack (OpenBlas) found?" FORCE)
        # set( lapack_FOUND TRUE CACHE BOOL "lapack (OpenBlas) found?" FORCE)
        set( lapack_flavor "openblas")
        set( lapack_url "http://www.openblas.net")
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

  if( NOT TARGET lapack AND NOT LAPACK_FOUND )
      MESSAGE( STATUS "Looking for lapack (no cmake config files)...")
      find_package( BLAS QUIET )

      if( BLAS_FOUND )
        find_package(LAPACK QUIET)
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

  set_package_properties( BLAS PROPERTIES
    URL "${lapack_url}"
    DESCRIPTION "Basic Linear Algebra Subprograms"
    TYPE OPTIONAL
    PURPOSE "Required for building the lapack_wrap component." )
  if( "${lapack_flavor}" STREQUAL "netlib")
    set_package_properties( lapack PROPERTIES
      URL "${lapack_url}"
      DESCRIPTION "Linear Algebra PACKage"
      TYPE OPTIONAL
      PURPOSE "Required for building the lapack_wrap component." )
  else()
    set_package_properties( LAPACK PROPERTIES
      URL "${lapack_url}"
      DESCRIPTION "Linear Algebra PACKage"
      TYPE OPTIONAL
      PURPOSE "Required for building the lapack_wrap component." )
  endif()
endmacro()

#------------------------------------------------------------------------------
# Helper macros for CUDA/Unix
#
# https://devblogs.nvidia.com/tag/cuda/
# https://devblogs.nvidia.com/building-cuda-applications-cmake/
#------------------------------------------------------------------------------
macro( setupCudaEnv )

  # if WITH_CUDA is set, use the provided value, otherwise disable CUDA unless
  # the CUDA_COMPILER exists.
  if( NOT DEFINED WITH_CUDA )
    if( EXISTS "${CMAKE_CUDA_COMPILER}" )
      set( WITH_CUDA ON)
    else()
      set( WITH_CUDA OFF)
    endif()
  endif()
  set( WITH_CUDA ${WITH_CUDA} CACHE BOOL "Attempt to compile CUDA kernels." )

  add_feature_info( Cuda WITH_CUDA "Build CUDA kernels for GPU compute.")

  if( WITH_CUDA AND NOT DEFINED CUDA_DBS_STRING )
    set( CUDA_DBS_STRING "CUDA" CACHE BOOL
      "If CUDA is available, this variable is 'CUDA'")

    set(OUTPUTFILE ${CMAKE_CURRENT_SOURCE_DIR}/config/cuda_script) # No suffix required
    set(CUDAFILE ${CMAKE_CURRENT_SOURCE_DIR}/config/query_gpu.cu)
    execute_process(COMMAND nvcc -lcuda ${CUDAFILE} -o ${OUTPUTFILE})
    execute_process(COMMAND ${OUTPUTFILE}
                    RESULT_VARIABLE CUDA_RETURN_CODE OUTPUT_VARIABLE ARCH)

    if (${CUDA_RETURN_CODE} EQUAL 0)
      message(STATUS "CUDA Architecture: ${ARCH}")
      # CMAKE currently only allows up to C++14 as the NVCC language level
      set(CMAKE_CUDA_STANDARD "14")
      set(CMAKE_CUDA_FLAGS "${ARCH} -g -G --expt-relaxed-constexpr" CACHE STRING "Standard CUDA flags"
        FORCE)
      set(CMAKE_CUDA_FLAGS_DEBUG "-O0" CACHE STRING "CUDA debug flags" FORCE)
      set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "-O2 --generate-line-info" CACHE
        STRING "CUDA release with debug information flags" FORCE)
      set(CMAKE_CUDA_FLAGS_RELEASE "-O2" CACHE STRING "CUDA release flags"
        FORCE)
    else()
      message(WARNING ${ARCH})
    endif()
  endif()

endmacro()

#------------------------------------------------------------------------------
# Setup QT (any)
#------------------------------------------------------------------------------
macro( setupQt )
  message( STATUS "Looking for Qt SDK...." )

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
    message( STATUS "Looking for Qt SDK....not found." )
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
      message( STATUS "Looking for Qt SDK....found ${QTDIR}." )
    else()
      set( QT_FOUND "QT-NOTFOUND" )
      message( STATUS "Looking for Qt SDK....not found." )
    endif()
  endif()

  if( QT_FOUND )
    mark_as_advanced( Qt5Core_DIR Qt5Gui_DIR Qt5Gui_EGL_LIBRARY
      Qt5Widgets_DIR QTDIR)
  endif()

  set_package_properties( Qt PROPERTIES
    URL "http://qt.io"
    DESCRIPTION "Qt is a comprehensive cross-platform C++ application framework."
    TYPE OPTIONAL
    PURPOSE "Only needed to demo qt version of draco_diagnostics." )

endmacro()

#------------------------------------------------------------------------------
# Setup GSL (any)
#------------------------------------------------------------------------------
macro( setupGSL )

  if( NOT TARGET GSL::gsl )

    message( STATUS "Looking for GSL..." )
    set( QUIET "QUIET")

    # There are 3 ways to find gsl:

    # 1. Config mode.
    #    If CMAKE_PREFIX_PATH contains a GSL install prefix directory and
    #    the file gsl-config.cmake is found somewhere in this installation
    #    tree, then the targets defined by gsl-config.cmake will be used.
    find_package( GSL CONFIG ${QUIET} )

  endif()

  if( NOT TARGET GSL::gsl ) # if option #1 was successful, skip this.

    # 2. pkg-config mode (Linux)
    #    IF GSL_ROOT_DIR isn't set, look for the binary 'gsl-config' in $PATH.
    #    If found, run it to discover and set GSL_ROOT_DIR that will be used
    #    in method #3.

    if( "$ENV{GSL_ROOT_DIR}x" STREQUAL "x" AND "${GSL_ROOT_DIR}x" STREQUAL "x")
      find_program( GSL_CONFIG gsl-config )
      if( EXISTS "${GSL_CONFIG}" )
        exec_program( "${GSL_CONFIG}"
          ARGS --prefix
          OUTPUT_VARIABLE GSL_ROOT_DIR )
      endif()
    endif()

    # 3. Module mode.
    #    Locate GSL by using the value of GSL_ROOT_DIR or by looking in
    #    standard locations. We add 'REQUIRED' here because if this fails,
    #    then we abort the built.
    find_package( GSL REQUIRED ${QUIET} )

  endif()

  # Print a report
  if( TARGET GSL::gsl )
    if( TARGET GSL::gsl AND NOT GSL_LIBRARY )
      foreach( config NOCONFIG DEBUG RELEASE RELWITHDEBINFO )
        get_target_property(tmp GSL::gsl IMPORTED_LOCATION_${config} )
        if( EXISTS ${tmp} AND NOT GSL_LIBRARY )
          set( GSL_LIBRARY ${tmp} )
        endif()
      endforeach()
    endif()
    message( STATUS "Looking for GSL.......found ${GSL_LIBRARY}" )
    mark_as_advanced( GSL_CONFIG_EXECUTABLE )
  else()
    message( STATUS "Looking for GSL.......not found" )
  endif()

  # If successful in finding GSL, provide some information for the vendor
  # summary reported by src/CMakeLists.txt.
  if( TARGET GSL::gsl )
    #=============================================================================
    # Include some information that can be printed by the build system.
    set_package_properties( GSL PROPERTIES
      URL "https://www.gnu.org/software/gsl"
      DESCRIPTION "The GNU Scientific Library (GSL) is a numerical library for C and C++
   programmers."
      TYPE REQUIRED
      PURPOSE "Required for rng and quadrature components." )
  endif()
  unset(QUIET)

endmacro()

#------------------------------------------------------------------------------
# Setup ParMETIS (any)
#------------------------------------------------------------------------------
macro( setupParMETIS )

  set( QUIET "QUIET")
  if( NOT TARGET METIS::metis )
    message( STATUS "Looking for METIS..." )

    find_package( METIS CONFIG ${QUIET} )
    if( NOT TARGET METIS::metis )
      find_package( METIS ${QUIET} )
    endif()
    if( TARGET METIS::metis )
      if( TARGET METIS::metis AND NOT METIS_LIBRARY )
        foreach( config NOCONFIG DEBUG RELEASE RELWITHDEBINFO )
          get_target_property(tmp METIS::metis IMPORTED_LOCATION_${config} )
          if( EXISTS ${tmp} AND NOT METIS_LIBRARY )
            set( METIS_LIBRARY ${tmp} )
          endif()
        endforeach()
      endif()
      message( STATUS "Looking for METIS.....found ${METIS_LIBRARY}" )
    else()
      message( STATUS "Looking for METIS.....not found" )
    endif()

    #=============================================================================
    # Include some information that can be printed by the build system.
    set_package_properties( METIS PROPERTIES
      DESCRIPTION "METIS"
      TYPE OPTIONAL
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
    else()
      message( STATUS "Looking for ParMETIS..not found" )
    endif()

    #=============================================================================
    # Include some information that can be printed by the build system.
    set_package_properties( ParMETIS PROPERTIES
      DESCRIPTION "MPI Parallel METIS"
      TYPE OPTIONAL
      URL "http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview"
      PURPOSE "ParMETIS is an MPI-based parallel library that implements a
   variety of algorithms for partitioning unstructured graphs, meshes, and for
   computing fill-reducing orderings of sparse matrices." )

  endif()
  unset(QUIET)
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
    else()
      message( STATUS "Looking for SuperLU_DIST.....not found" )
    endif()

    if( ${SuperLU_DIST_VERSION} VERSION_GREATER 5.2.9 )
      message( FATAL_ERROR "The API change in SuperLU_DIST 5.3+ is not yet
      supported by Draco. Please use a version of SuperLU_DIST prior to 5.3.")
    endif()

    #===========================================================================
    # Include some information that can be printed by the build system.
    set_package_properties( SuperLU_DIST PROPERTIES
      URL " http://crd-legacy.lbl.gov/~xiaoye/SuperLU/"
      DESCRIPTION "SuperLU_DIST"
      TYPE OPTIONAL
      PURPOSE "SuperLU is a general purpose library for the direct solution of
   large, sparse, nonsymmetric systems of linear equations on high performance
   machines."  )

  endif()

endmacro()

#------------------------------------------------------------------------------
# Setup Eospac (https://laws.lanl.gov/projects/data/eos.html)
#------------------------------------------------------------------------------
macro( setupEOSPAC )

  if( NOT TARGET EOSPAC::eospac )
    message( STATUS "Looking for EOSPAC..." )

    find_package( EOSPAC QUIET )

    if( EOSPAC_FOUND )
      message( STATUS "Looking for EOSPAC....found ${EOSPAC_LIBRARY}" )
    else()
      message( STATUS "Looking for EOSPAC....not found" )
    endif()

    #===========================================================================
    # Include some information that can be printed by the build system.
    set_package_properties( EOSPAC PROPERTIES
      URL "https://laws.lanl.gov/projects/data/eos.html"
      DESCRIPTION "Access SESAME thermodynamic and transport data."
      TYPE OPTIONAL
      PURPOSE "Required for bulding the cdi_eospac component." )
  endif()

endmacro()

#------------------------------------------------------------------------------
# Setup COMPTON (https://gitlab.lanl.gov/keadyk/CSK_generator)
#------------------------------------------------------------------------------
macro( setupCOMPTON )

  if( NOT TARGET compton::compton )
    message( STATUS "Looking for COMPTON..." )

    find_package( COMPTON QUIET )

    if( COMPTON_FOUND )
      message( STATUS "Looking for COMPTON...found ${COMPTON_LIBRARY}" )
    else()
      message( STATUS "Looking for COMPTON...not found" )
    endif()

    #===========================================================================
    # Include some information that can be printed by the build system.
    set_package_properties( COMPTON PROPERTIES
      URL "https://gitlab.lanl.gov/CSK/CSK"
      DESCRIPTION "Access multigroup Compton scattering data."
      TYPE OPTIONAL
      PURPOSE "Required for bulding the compton component." )
  endif()

endmacro()

#------------------------------------------------------------------------------
# Helper macros for setup_global_libraries()
#------------------------------------------------------------------------------
macro( SetupVendorLibrariesUnix )

  setupGSL()
  setupParMETIS()
  setupSuperLU_DIST()
  setupCOMPTON()
  setupEospac()
  setupRandom123()
  setupCudaEnv()
  setupPython()
  setupQt()

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

  # Doxygen ------------------------------------------------------------------
  message( STATUS "Looking for Doxygen..." )
  find_package( Doxygen QUIET OPTIONAL_COMPONENTS dot mscgen dia )
  set_package_properties( Doxygen PROPERTIES
    URL "http://www.stack.nl/~dimitri/doxygen"
    DESCRIPTION "Doxygen autodoc generator"
    TYPE OPTIONAL
    PURPOSE "Required for building develop HTML documentation."
    )
  if( DOXYGEN_FOUND )
    message( STATUS "Looking for Doxygen...found version ${DOXYGEN_VERSION}" )
  else()
    message( STATUS "Looking for Doxygen...not found" )
  endif()

endmacro()

##---------------------------------------------------------------------------##
## Vendors for building on Windows-based platforms.
##---------------------------------------------------------------------------##

macro( SetupVendorLibrariesWindows )

  setupGSL()
  setupParMETIS()
  setupSuperLU_DIST()
  setupRandom123()
  setupCOMPTON()
  setupEospac()
  setupPython()
  setupQt()
  setupCudaEnv()

  # Doxygen ------------------------------------------------------------------
  message( STATUS "Looking for Doxygen..." )
  find_package( Doxygen QUIET OPTIONAL_COMPONENTS dot mscgen dia )
  set_package_properties( Doxygen PROPERTIES
    URL "http://www.stack.nl/~dimitri/doxygen"
    DESCRIPTION "Doxygen autodoc generator"
    TYPE OPTIONAL
    PURPOSE "Required for building develop HTML documentation."
    )
  if( DOXYGEN_FOUND )
    message( STATUS "Looking for Doxygen...found version ${DOXYGEN_VERSION}" )
  else()
    message( STATUS "Looking for Doxygen...not found" )
  endif()


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
  if( NOT DEFINED VENDOR_DIR AND IS_DIRECTORY "$ENV{VENDOR_DIR}" )
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
  endif()

  # Import environment variables related to vendors
  # 1. Use command line variables (-DLAPACK_LIB_DIR=<path>
  # 2. Use environment variables ($ENV{LAPACK_LIB_DIR}=<path>)
  # 3. Try to find vendor in $VENDOR_DIR
  # 4. Don't set anything and let the user set a value in the cache
  #    after failed 1st configure attempt.
  if( NOT DEFINED LAPACK_LIB_DIR AND IS_DIRECTORY $ENV{LAPACK_LIB_DIR} )
    set( LAPACK_LIB_DIR $ENV{LAPACK_LIB_DIR} )
    set( LAPACK_INC_DIR $ENV{LAPACK_INC_DIR} )
  endif()
  if( NOT LAPACK_LIB_DIR AND IS_DIRECTORY ${VENDOR_DIR}/lapack-3.4.2/lib )
    set( LAPACK_LIB_DIR "${VENDOR_DIR}/lapack-3.4.2/lib" )
    set( LAPACK_INC_DIR "${VENDOR_DIR}/lapack-3.4.2/include" )
  endif()

  if( NOT DEFINED GSL_LIB_DIR )
    if( IS_DIRECTORY $ENV{GSL_LIB_DIR}  )
      set( GSL_LIB_DIR $ENV{GSL_LIB_DIR} )
      set( GSL_INC_DIR $ENV{GSL_INC_DIR} )
    elseif( IS_DIRECTORY ${VENDOR_DIR}/gsl/lib )
      set( GSL_LIB_DIR "${VENDOR_DIR}/gsl/lib" )
      set( GSL_INC_DIR "${VENDOR_DIR}/gsl/include" )
    endif()
  endif()

  if( NOT DEFINED ParMETIS_ROOT_DIR )
    if( IS_DIRECTORY $ENV{ParMETIS_ROOT_DIR}  )
      set( ParMETIS_ROOT_DIR $ENV{ParMETIS_ROOT_DIR} )
    endif()
  endif()

  if( NOT DEFINED RANDOM123_INC_DIR AND IS_DIRECTORY $ENV{RANDOM123_INC_DIR}  )
    set( RANDOM123_INC_DIR $ENV{RANDOM123_INC_DIR} )
  endif()
  if( NOT DEFINED RANDOM123_INC_DIR AND
      IS_DIRECTORY ${VENDOR_DIR}/Random123-1.08/include )
    set( RANDOM123_INC_DIR "${VENDOR_DIR}/Random123-1.08/include" )
  endif()

endmacro()

#------------------------------------------------------------------------------
# This macro should contain all the system libraries which are required to link
# the main objects.
# ------------------------------------------------------------------------------
macro( setupVendorLibraries )

  message( "\nVendor Setup:\n")

  #
  # General settings
  #
  setVendorVersionDefaults()
  if( NOT TARGET lapack )
    setupLAPACKLibraries()
  endif()

  # System specific settings
  if ( UNIX )
    setupMPILibrariesUnix()
    setupVendorLibrariesUnix()
  elseif( WIN32 )
    setupMPILibrariesWindows()
    setupVendorLibrariesWindows()
  else()
    message( FATAL_ERROR "
I don't know how to setup global (vendor) libraries for this platform.
WIN32=0; UNIX=0; CMAKE_SYSTEM=${CMAKE_SYSTEM};
CMAKE_SYSTEM_NAME=${CMAKE_SYSTEM_NAME}" )
  endif()

  # Add commands to draco-config.cmake (which is installed for use by othe
  # projects), to setup Draco's vendors
  set( Draco_EXPORT_TARGET_PROPERTIES "${Draco_EXPORT_TARGET_PROPERTIES}
macro( dbs_basic_setup )

  message(\"
Looking for Draco...\")
  message(\"Looking for Draco...\${draco_DIR}
  \")

  # Provide helper functions used by component CMakeLists.txt files
  # This block of code generated by draco/config/vendor_libraries.cmake.

  # CMake macros that check the system for features like 'gethostname', etc.
  include( platform_checks )

  # Sanity check for Cray Programming Environments
  query_craype()

  # Set compiler options
  include( compilerEnv )
  dbsSetupCxx()
  dbsSetupFortran()
  dbsSetupProfilerTools()

  # CMake macros like 'add_component_library' and 'add_component_executable'
  include( component_macros )

  # CMake macros to query the availability of TPLs.
  include( vendor_libraries )

  # Provide targets for MPI, Metis, etc.
  setupVendorLibraries()

endmacro()
")

  message( " " )

endmacro()

#----------------------------------------------------------------------#
# End vendor_libraries.cmake
#----------------------------------------------------------------------#
