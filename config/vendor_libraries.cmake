#-----------------------------*-cmake-*----------------------------------------#
# file   config/global_libraries.cmake
# author 
# date   2010 June 6
# brief  Setup Vendors
# note   © Copyright 2010 LANS, LLC  
#------------------------------------------------------------------------------#
# $Id$ 
#------------------------------------------------------------------------------#

#
# Look for any libraries which are required at the toplevel.
# 

#------------------------------------------------------------------------------
# Helper macros for setup_global_libraries()
# Assign here the library version to be used.
#------------------------------------------------------------------------------
macro( setVendorVersionDefaults )

  # set( TSP_BOOST_VERSION 1.40.0 CACHE STRING
    # "Minimum supported Boost version." )
  # set( XERCESC_VERSION 2.8.0 )
  # if( NOT ONLY_LIBTRANSPIRE )
    # set( LAPACK_VERSION 3.1.1 ) 
    # set( HDF5_VERSION 1.8.0 )
    # set( ALM_VERSION "" )
    # set( EXPAT_VERSION "" )
    # set( OPENSSL_VERSION "" )
    # set( SIMMETRIX_VERSION "6.3-080514" )
    # set( BISON_VERSION "" )
    # set( FLEX_VERSION "" )
  # endif( NOT ONLY_LIBTRANSPIRE )

  #Set the preferred search directories(ROOT)

  #Check that VENDOR_DIR is defined as a cache variable or as an
  #environment variable. If defined as both then take the
  #environment variable.

  # See if VENDOR_DIR is set.  Try some defaults if it is not set.
  if( NOT VENDOR_DIR AND IS_DIRECTORY "$ENV{VENDOR_DIR}" )
    set( VENDOR_DIR $ENV{VENDOR_DIR} )
  endif()
  # If needed, try some obvious palces.
  if( NOT VENDOR_DIR )
     if( IS_DIRECTORY /ccs/codes/radtran/vendors/Linux64 )
        set( VENDOR_DIR /ccs/codes/radtran/vendors/Linux64 )
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
    set( VENDOR_DIR $ENV{VENDOR_DIR} CACHE PATH
      "Root directory where Transpire 3rd party libraries are located." )
  else( IS_DIRECTORY "${VENDOR_DIR}")
    message( "
WARNING: VENDOR_DIR not defined locally or in user environment,
individual vendor directories should be defined." )
  endif( IS_DIRECTORY "${VENDOR_DIR}")

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
  if( NOT LAPACK_LIB_DIR AND IS_DIRECTORY ${VENDOR_DIR}/clapack/lib )
     set( LAPACK_LIB_DIR "${VENDOR_DIR}/clapack/lib" )
     set( LAPACK_INC_DIR "${VENDOR_DIR}/clapack/include" )
  endif()

  if( NOT GSL_LIB_DIR AND IS_DIRECTORY $ENV{GSL_LIB_DIR}  )
     set( GSL_LIB_DIR $ENV{GSL_LIB_DIR} )
     set( GSL_INC_DIR $ENV{GSL_INC_DIR} )
  endif()
  if( NOT GSL_LIB_DIR AND IS_DIRECTORY ${VENDOR_DIR}/gsl/lib )
     set( GSL_LIB_DIR "${VENDOR_DIR}/gsl/lib" )
     set( GSL_INC_DIR "${VENDOR_DIR}/gsl/include" )
  endif()

  #Set the preferred search paths
  
  # Look for special build of boost (/DSECURE_SCL=0).
  # set( BOOST_ROOT ${VENDOR_DIR}/boost-${TSP_BOOST_VERSION}-sscl0 )
  # if( NOT EXISTS ${BOOST_ROOT} )
    # set( BOOST_ROOT ${VENDOR_DIR}/boost-${TSP_BOOST_VERSION} )
  # endif()
  # set( BOOST_ROOT ${BOOST_ROOT} 
      # CACHE PATH "Where should CMake loook for Boost?" )
      
  # set( XERCESC_ROOT
    # ${VENDOR_DIR}/xerces-c-${XERCESC_VERSION} )
        
    # if    (NOT ONLY_LIBTRANSPIRE )
      # set( LAPACK_ROOT  ${VENDOR_DIR}/lapack-${LAPACK_VERSION} )
      # set( HDF5_ROOT    ${VENDOR_DIR}/hdf5-${HDF5_VERSION} )
      # set( ALM_ROOT 
        # ${transpire_SOURCE_DIR}/lib/alm/${CMAKE_SYSTEM_PROCESSOR}-${CMAKE_SYSTEM_NAME} )
      # set( EXPAT_ROOT   ${VENDOR_DIR}/Expat )
      # set( OPENSSL_ROOT ${VENDOR_DIR}/OpenSSL )
      # set( SIMMETRIX_ROOT ${VENDOR_DIR}/simmetrix-${SIMMETRIX_VERSION} )
      # set( BISON_ROOT   ${VENDOR_DIR}/bison-${BISON_VERSION} )
      # set( FLEX_ROOT    ${VENDOR_DIR}/flex-${FLEX_VERSION} )
    # endif (NOT ONLY_LIBTRANSPIRE )

    # if   ( VERBOSE )
      # message("")
      # message("Preferred search directories:")
      # message("BOOST_ROOT     = ${BOOST_ROOT}") 
      # message("XERCESC_ROOT   = ${XERCESC_ROOT}")
      # if    (NOT ONLY_LIBTRANSPIRE )
        # message("LAPACK_ROOT    = ${LAPACK_ROOT}") 
        # message("HDF5_ROOT      = ${HDF5_ROOT}") 
        # message("ALM_ROOT       = ${ALM_ROOT}") 
        # message("EXPAT_ROOT     = ${EXPAT_ROOT}") 
        # message("OPENSSL_ROOT   = ${OPENSSL_ROOT}") 
        # message("SIMMETRIX_ROOT = ${SIMMETRIX_ROOT}") 
        # message("BISON_ROOT     = ${BISON_ROOT}") 
        # message("FLEX_ROOT      = ${FLEX_ROOT}") 
      # endif (NOT ONLY_LIBTRANSPIRE )
      # message("--VERBOSE--")
      # message("")
    # endif( VERBOSE )

endmacro()

#------------------------------------------------------------------------------
# Helper macros for setup_global_libraries()
#------------------------------------------------------------------------------
macro( SetupVendorLibrariesWindows )

# This module will set the following variables:
#   MPI_FOUND                  TRUE if we have found MPI
#   MPI_COMPILE_FLAGS          Compilation flags for MPI programs
#   MPI_INCLUDE_PATH           Include path(s) for MPI header
#   MPI_LINK_FLAGS             Linking flags for MPI programs
#   MPI_LIBRARY                First MPI library to link against (cached)
#   MPI_EXTRA_LIBRARY          Extra MPI libraries to link against (cached)
#   MPI_LIBRARIES              All libraries to link MPI programs against
#   MPIEXEC                    Executable for running MPI programs
#   MPIEXEC_NUMPROC_FLAG       Flag to pass to MPIEXEC before giving it the
#                              number of processors to run on
#   MPIEXEC_PREFLAGS           Flags to pass to MPIEXEC directly before the
#                              executable to run.
#   MPIEXEC_POSTFLAGS          Flags to pass to MPIEXEC after all other flags.

  if( IS_DIRECTORY "$ENV{MPI_INC_DIR}" AND IS_DIRECTORY "$ENV{MPI_LIB_DIR}" )
     set( MPI_INCLUDE_PATH $ENV{MPI_INC_DIR} )
     set( MPI_LIBRARY $ENV{MPI_LIB_DIR}/mpi.lib )
  endif()
  find_package( MPI )
  if( MPI_FOUND )
     set( DRACO_C4 "MPI" )  
  else()
     set( DRACO_C4 "SCALAR" )
  endif()
  set( DRACO_C4 "${DRACO_C4}" CACHE STRING "C4 communication mode (SCALAR or MPI)" )
  if( "${DRACO_C4}" STREQUAL "MPI"    OR 
      "${DRACO_C4}" STREQUAL "SCALAR" )
          # message("
          # MPI_FOUND         = ${MPI_FOUND}
          # MPI_LIBRARY       = ${MPI_LIBRARY}
          # MPI_EXTRA_LIBRARY = ${MPI_EXTRA_LIBRARY}
          # MPI_LIBRARIES     = ${MPI_LIBRARIES}
          # MPIEXEC           = ${MPIEXEC}
          # MPIEXEC_NUMPROC_FLAG = ${MPIEXEC_NUMPROC_FLAG}
          # DRACO_C4          = ${DRACO_C4}
          # ")
  else()
    message( FATAL_ERROR "DRACO_C4 must be either MPI or SCALAR" )
  endif()
    
  # LAPACK & BLAS
  find_package( CLAPACK REQUIRED )
  
  # GSL
  find_package( GSL REQUIRED )
  
  
  # find_package( XercesC QUIET )

  # Boost: 
  #
  # required by the build system (not by any files that are
  # released). 
  # set( Boost_USE_STATIC_LIBS ON CACHE BOOL 
       # "If available use .lib libs instead of .dll.")
  #set( Boost_DEBUG ON )
  # find_package( Boost ${TSP_BOOST_VERSION}
    # REQUIRED
      # date_time 
      # filesystem 
      # program_options
      # regex  
      # system
      # test_exec_monitor
      # unit_test_framework
    # )
  # add_definitions( ${Boost_LIB_DIAGNOSTIC_DEFINITIONS} )

# Haven't bothered to try OpenCASCADE on w32 so we just set
# this option to OFF.
  # option( ENABLE_OCC "Enable programs linked against OpenCASCADE" OFF )

  # Use Qt 4 if we want.
  # option( ENABLE_QT4 "Enable program linked against Qt 4" OFF )

  # if ( ENABLE_QT4 )
    # find_package( Qt4 )
  # endif( ENABLE_QT4 )

endmacro()

#------------------------------------------------------------------------------
# Helper macros for setup_global_libraries()
#------------------------------------------------------------------------------
macro( SetupVendorLibrariesUnix )

   # This module will set the following variables:
   #   MPI_FOUND                  TRUE if we have found MPI
   #   MPI_COMPILE_FLAGS          Compilation flags for MPI programs
   #   MPI_INCLUDE_PATH           Include path(s) for MPI header
   #   MPI_LINK_FLAGS             Linking flags for MPI programs
   #   MPI_LIBRARY                First MPI library to link against (cached)
   #   MPI_EXTRA_LIBRARY          Extra MPI libraries to link against (cached)
   #   MPI_LIBRARIES              All libraries to link MPI programs against
   #   MPIEXEC                    Executable for running MPI programs
   #   MPIEXEC_NUMPROC_FLAG       Flag to pass to MPIEXEC before giving it the
   #                              number of processors to run on
   #   MPIEXEC_PREFLAGS           Flags to pass to MPIEXEC directly before the
   #                              executable to run.
   #   MPIEXEC_POSTFLAGS          Flags to pass to MPIEXEC after all other flags.
   
   if( IS_DIRECTORY "$ENV{MPI_INC_DIR}" AND IS_DIRECTORY "$ENV{MPI_LIB_DIR}" )
      set( MPI_INCLUDE_PATH $ENV{MPI_INC_DIR} )
      set( MPI_LIBRARY $ENV{MPI_LIB_DIR}/libmpi.so )
      set( MPI_EXTRA_LIBRARY $ENV{MPI_LIB_DIR}/libmpi_cxx.so )
   endif()
   find_package( MPI )
   if( MPI_FOUND )
      set( DRACO_C4 "MPI" )  
      set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DOMPI_SKIP_MPICXX" )
   else()
      set( DRACO_C4 "SCALAR" )
   endif()
   set( DRACO_C4 "${DRACO_C4}" CACHE STRING "C4 communication mode (SCALAR or MPI)" )
   if( "${DRACO_C4}" STREQUAL "MPI"    OR 
       "${DRACO_C4}" STREQUAL "SCALAR" )
      # message("
      #     MPI_FOUND         = ${MPI_FOUND}
      #     MPI_LIBRARY       = ${MPI_LIBRARY}
      #     MPI_EXTRA_LIBRARY = ${MPI_EXTRA_LIBRARY}
      #     MPI_LIBRARIES     = ${MPI_LIBRARIES}
      #     MPIEXEC           = ${MPIEXEC}
      #     MPIEXEC_NUMPROC_FLAG = ${MPIEXEC_NUMPROC_FLAG}
      #     DRACO_C4          = ${DRACO_C4}
      #     ")
   else()
      message( FATAL_ERROR "DRACO_C4 must be either MPI or SCALAR" )
   endif()

  # LAPACK & BLAS
  find_package( CLAPACK REQUIRED )
 
  # GSL
  find_package( GSL REQUIRED )

  # Gandolf
  find_package( Gandolf REQUIRED )

  # find_package( XercesC REQUIRED )

  # Boost: 
  #
  # required by the build system (not by any files that are
  # released). See libtranspire/generate_build_date_f90.cpp.

  # set( Boost_USE_MULTITHREADED ON )
  # set( Boost_USE_STATIC_LIBS ON CACHE BOOL 
       # "If available use .lib libs instead of .so.")

#  set( Boost_DEBUG ON )
  # find_package( Boost ${TSP_BOOST_VERSION}
    # REQUIRED
      # date_time
      # filesystem
      # program_options
      # regex
      # system 
      # test_exec_monitor
      # unit_test_framework
    # )

  # If the user has selected to build the OpenCASCADE programs, then
  # try to track down the OpenCASCADE installation.
  # option( ENABLE_OCC "Enable programs linked against OpenCASCADE" OFF )

  # if ( ENABLE_OCC )
    # find_package( OCC )
    # find_package( SQLITE3 )
  # endif( ENABLE_OCC )

  # Use Qt 4 if we want.
  # option( ENABLE_QT4 "Enable program linked against Qt 4" OFF )

  # if ( ENABLE_QT4 )
    # find_package( Qt4 )
  # endif( ENABLE_QT4 )

endmacro()

#------------------------------------------------------------------------------
# This macro should contain all the system libraries which are
# required to link the main objects.
#------------------------------------------------------------------------------
macro( setupVendorLibraries )

  #
  # General settings
  # 
  setVendorVersionDefaults()

  # System specific settings
  if ( UNIX )
    setupVendorLibrariesUnix()
  elseif( WIN32 )
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
