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

#------------------------------------------------------------------------------#
# Setup MPI when on Linux
#------------------------------------------------------------------------------#
macro( setupMPILibrariesUnix )

   # MPI ---------------------------------------------------------------------
   if( NOT "${DRACO_C4}" STREQUAL "SCALAR" )
      # Try to find MPI in the default locations (look for mpic++ in PATH)
      # This module will set the following variables:
      #   MPI_FOUND                  TRUE have we found MPI
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
      
      message(STATUS "Looking for MPI...")

      # Preserve data that may already be set.
      if( DEFINED MPIEXEC )
         set( MPIEXEC ${MPIEXEC} CACHE FILEPATH "Executable for running MPI programs." )
      endif()
      if( DEFINED MPIEXEC_MAX_NUMPROCS )
         set( MPIEXEC_MAX_NUMPROCS ${MPIEXEC_MAX_NUMPROCS} CACHE STRING
            "Maximum number of processors available to run MPI applications.")
      endif()
      if( DEFINED MPIEXEC_NUMPROC_FLAG )
         set( MPIEXEC_NUMPROC_FLAG ${MPIEXEC_NUMPROC_FLAG} CACHE
            STRING "Flag used by MPI to specify the number of processes for MPIEXEC; the next option will be the number of processes." )
      endif()

      # Before calling find_package(MPI), if CMAKE_<lang>_COMPILER is
      # actually an MPI wrapper script:
      #
      # - set MPI_<lang>_COMPILER, so that FindMPI will skip automatic
      #   compiler detection; and
      # - set MPI_<lang>_NO_INTERROGATE, so that FindMPI will skip MPI
      #   include-path and library detection and trust that the
      #   provided compilers can build MPI binaries without further
      #   help.
      get_filename_component( compiler_wo_path
         "${CMAKE_CXX_COMPILER}" NAME )
      if( "${compiler_wo_path}" MATCHES "mpi" )
         set( MPI_CXX_COMPILER ${CMAKE_CXX_COMPILER} )
         set( MPI_CXX_NO_INTERROGATE ${CMAKE_CXX_COMPILER} )
      endif()

      get_filename_component( compiler_wo_path
         "${CMAKE_C_COMPILER}" NAME )
      if( "${compiler_wo_path}" MATCHES "mpi" )
         set( MPI_C_COMPILER ${CMAKE_C_COMPILER} )
         set( MPI_C_NO_INTERROGATE ${CMAKE_C_COMPILER} )
      endif()

      get_filename_component( compiler_wo_path
         "${CMAKE_Fortran_COMPILER}" NAME )
      if( "${compiler_wo_path}" MATCHES "mpi" )
         set( MPI_Fortran_COMPILER ${CMAKE_Fortran_COMPILER} )
         set( MPI_Fortran_NO_INTERROGATE ${CMAKE_Fortran_COMPILER} )
      endif()

      # On CT, only look for .a MPI libraries when build type is STATIC
      # if( "${SITE}" MATCHES "c[it]" AND
      #       "${DRACO_LIBRARY_TYPE}" MATCHES "STATIC" )
      #    set( CMAKE_FIND_LIBRARY_SUFFIXES .a ) # automatic for Catamount.cmake
      # endif()

      # First attempt to find mpi
      find_package( MPI QUIET )

      # On CI/CT, we need to link to STATIC libraries, but
      # findMPI.cmake will conly return shared libraries.  Force
      # conversion to static iff CT with STATIC and openmpi are
      # selected.
#       if( "${SITE}" MATCHES "c[it]" AND
#             "${MPIEXEC}" MATCHES "openmpi" AND
#             "${DRACO_LIBRARY_TYPE}" MATCHES "STATIC" )
#          foreach( VAR MPI_CXX_LIBRARIES MPI_C_LIBRARIES 
#                MPI_Fortran_LIBRARIES MPI_EXTRA_LIBRARY MPI_LIBRARIES )
#             string( REGEX REPLACE "mpi_cxx[.]so" "mpi_cxx.a" ${VAR} "${${VAR}}" )
#             string( REGEX REPLACE "mpi[.]so" "mpi.a" ${VAR} "${${VAR}}" )
#             string( REGEX REPLACE "mpi_mpifh[.]so" "mpi_mpifh.a" ${VAR} "${${VAR}}" )
#             string( REGEX REPLACE "mpi_usempif08[.]so" "mpi_usempif08.a" ${VAR} "${${VAR}}" )
#             string( REGEX REPLACE "pmi[.]so" "pmi.a" ${VAR} "${${VAR}}" )
#             set( ${VAR} "${${VAR}}" CACHE STRING
#                "MPI libraries to link against." FORCE )
#       endif()

      # Set Draco build system variables based on what we know about MPI.
      if( MPI_FOUND )
         set( DRACO_C4 "MPI" )  
         if( NOT MPIEXEC )
            if( "${SITE}" MATCHES "c[it]" ) # AND "${MPIEXEC}" MATCHES "openmpi"
               set( MPIEXEC aprun )
            else()
               message( FATAL_ERROR 
                  "MPI found but mpirun not in PATH. Aborting" )
            endif()
         endif()
         # Try to find the fortran mpi library
         if( EXISTS ${MPI_LIB_DIR} )
            find_library( MPI_Fortran_LIB mpi_f77 HINTS ${MPI_LIB_DIR} )
            mark_as_advanced( MPI_Fortran_LIB )
         endif()
         set( MPI_LIBRARIES "${MPI_LIBRARIES}" CACHE FILEPATH 
            "no mpi library for scalar build." FORCE )
         set( MPI_Fortran_LIBRARIES "${MPI_Fortran_LIBRARIES}" CACHE FILEPATH 
            "no mpi library for scalar build." FORCE )
         
      else()
         set( DRACO_C4 "SCALAR" )
         set( MPI_INCLUDE_PATH "" CACHE FILEPATH 
            "no mpi library for scalar build." FORCE )
         set( MPI_LIBRARY "" CACHE FILEPATH 
            "no mpi library for scalar build." FORCE )
         set( MPI_LIBRARIES "" CACHE FILEPATH 
            "no mpi library for scalar build." FORCE )
      endif()

      # Save the result in the cache file.
      set( DRACO_C4 "${DRACO_C4}" CACHE STRING 
         "C4 communication mode (SCALAR or MPI)" )
      # Provide a constrained pull down list in cmake-gui
      set_property( CACHE DRACO_C4 PROPERTY STRINGS SCALAR MPI )
      if( "${DRACO_C4}" STREQUAL "MPI"    OR 
            "${DRACO_C4}" STREQUAL "SCALAR" )
         # no-op
      else()
         message( FATAL_ERROR "DRACO_C4 must be either MPI or SCALAR" )
      endif()

      execute_process( COMMAND ${MPIEXEC} --version
         OUTPUT_VARIABLE DBS_MPI_VER_OUT 
         ERROR_VARIABLE DBS_MPI_VER_ERR)

      set( DBS_MPI_VER "${DBS_MPI_VER_OUT} ${DBS_MPI_VER_ERR}") 

      # Check flavor and add optional flags
      if( "${MPIEXEC}" MATCHES openmpi OR
            "${DBS_MPI_VER}" MATCHES open-mpi )
         set( MPI_FLAVOR "openmpi" CACHE STRING "Flavor of MPI." )

         # Find the version of OpenMPI

         if( "${DBS_MPI_VER}" MATCHES "[0-9].[0-9].[0-9]" )
            string( REGEX REPLACE ".*([0-9]).([0-9]).([0-9]).*" "\\1"
               DBS_MPI_VER_MAJOR ${DBS_MPI_VER} )
            string( REGEX REPLACE ".*([0-9]).([0-9]).([0-9]).*" "\\2"
               DBS_MPI_VER_MINOR ${DBS_MPI_VER} )
         elseif( "${DBS_MPI_VER}" MATCHES "[0-9].[0-9]" )
            string( REGEX REPLACE ".*([0-9]).([0-9]).*" "\\1"
               DBS_MPI_VER_MAJOR ${DBS_MPI_VER} )
            string( REGEX REPLACE ".*([0-9]).([0-9]).*" "\\2"
               DBS_MPI_VER_MINOR ${DBS_MPI_VER} )
         endif()

         # sanity check, these OpenMPI flags (below) require version >= 1.4
         if( ${DBS_MPI_VER_MAJOR}.${DBS_MPI_VER_MINOR} VERSION_LESS 1.4 )
            message( FATAL_ERROR "OpenMPI version < 1.4 found." )
         endif()
               
         # Find cores/cpu and cpu/node.

         set( MPI_CORES_PER_CPU 4 )
         set( MPI_PHYSICAL_CORES 0 )
         if( EXISTS "/proc/cpuinfo" )
            file( READ "/proc/cpuinfo" cpuinfo_data )
            string( REGEX REPLACE "\n" ";" cpuinfo_data "${cpuinfo_data}" )
            foreach( line ${cpuinfo_data} )
               if( "${line}" MATCHES "cpu cores" )
                  string( REGEX REPLACE ".* ([0-9]+).*" "\\1"
                     MPI_CORES_PER_CPU "${line}" )
               elseif( "${line}" MATCHES "physical id" )
                  string( REGEX REPLACE ".* ([0-9]+).*" "\\1" tmp "${line}" )
                  if( ${tmp} GREATER ${MPI_PHYSICAL_CORES} )
                     set( MPI_PHYSICAL_CORES ${tmp} )
                  endif()
               endif()
            endforeach()
         endif()
         math( EXPR MPI_CPUS_PER_NODE "${MPIEXEC_MAX_NUMPROCS} / ${MPI_CORES_PER_CPU}" )
         set( MPI_CPUS_PER_NODE ${MPI_CPUS_PER_NODE} CACHE STRING
            "Number of multi-core CPUs per node" FORCE )
         set( MPI_CORES_PER_CPU ${MPI_CORES_PER_CPU} CACHE STRING
            "Number of cores per cpu" FORCE )

         # Check for hyperthreading - This is important for reserving
         # threads for OpenMP tests...

         # correct base-zero indexing
         math( EXPR MPI_PHYSICAL_CORES "${MPI_PHYSICAL_CORES} + 1" )
         math( EXPR MPI_MAX_NUMPROCS_PHYSICAL
            "${MPI_PHYSICAL_CORES} * ${MPI_CORES_PER_CPU}" )
         if( "${MPI_MAX_NUMPROCS_PHYSICAL}" STREQUAL "${MPIEXEC_MAX_NUMPROCS}" )
            set( MPI_HYPERTHREADING "OFF" CACHE BOOL "Are we using hyperthreading?" FORCE )
         else()
            set( MPI_HYPERTHREADING "ON" CACHE BOOL "Are we using hyperthreading?" FORCE )
         endif()
         
         #
         # EAP's flags can be found in Test.rh/General/run_job.pl
         # (look for $other_args).  In particular, it may be useful to
         # examine EAP's options for srun or aprun.
         # 
         # \sa
         # http://blogs.cisco.com/performance/open-mpi-v1-5-processor-affinity-options/
         #
         # --cpus-per-proc <N> option is needed for multi-threaded
         #   processes.  Use this as "threads per MPI rank."
         # --bind-to-socket will bind MPI ranks ordered by socket first (rank 0 onto
         #   socket 0, then rank 1 onto socket 1, etc.)

         # --bind-to-core added in OpenMPI-1.4
         set( MPIEXEC_OMP_POSTFLAGS 
            -bind-to-socket -cpus-per-proc ${MPI_CORES_PER_CPU} --report-bindings
            CACHE STRING "extra mpirun flags (list)." FORCE )
         set( MPIEXEC_OMP_POSTFLAGS_STRING 
            "-bind-to-socket -cpus-per-proc ${MPI_CORES_PER_CPU} --report-bindings"
            CACHE STRING "extra mpirun flags (list)." FORCE)
         mark_as_advanced( MPI_FLAVOR MPIEXEC_OMP_POSTFLAGS_STRING MPIEXEC_OMP_POSTFLAGS
            MPI_LIBRARIES )         

      elseif( "${MPIEXEC}" MATCHES aprun)
         # According to email from Mike McKay (2013/04/18), we might
         # need to set the the mpirun command to something like: 
         #
         #   setenv OMP_NUM_THREADS ${MPI_CORES_PER_CPU}; aprun ... -d ${MPI_CORES_PER_CPU}
         # 
         # consider '-cc none'
         if( NOT "$ENV{OMP_NUM_THREADS}x" STREQUAL "x" )
            set( MPI_CPUS_PER_NODE 1 CACHE STRING
               "Number of multi-core CPUs per node" FORCE )
            set( MPI_CORES_PER_CPU $ENV{OMP_NUM_THREADS} CACHE STRING
               "Number of cores per cpu" FORCE )
            set( MPIEXEC_OMP_POSTFLAGS -d ${MPI_CORES_PER_CPU} CACHE
               STRING "extra mpirun flags (list)." FORCE)
            set( MPIEXEC_OMP_POSTFLAGS_STRING "-d ${MPI_CORES_PER_CPU}" CACHE
               STRING "extra mpirun flags (string)." FORCE)
         else()
            message( STATUS "
WARNING: ENV{OMP_NUM_THREADS} is not set in your environment, 
         all OMP tests will be disabled." )
         endif()

      elseif( "${MPIEXEC}" MATCHES srun)
         if( NOT "$ENV{OMP_NUM_THREADS}x" STREQUAL "x" )
            set( MPI_CPUS_PER_NODE 1 CACHE STRING
               "Number of multi-core CPUs per node" FORCE )
            set( MPI_CORES_PER_CPU $ENV{OMP_NUM_THREADS} CACHE STRING
               "Number of cores per cpu" FORCE )
            set( MPIEXEC_OMP_POSTFLAGS -c${MPI_CORES_PER_CPU} CACHE
               STRING "extra mpirun flags (list)." FORCE)
            set( MPIEXEC_OMP_POSTFLAGS_STRING "-c${MPI_CORES_PER_CPU}" CACHE
               STRING "extra mpirun flags (string)." FORCE)
         else()
            message( STATUS "
WARNING: ENV{OMP_NUM_THREADS} is not set in your environment, 
         all OMP tests will be disabled." )
         endif()

         set( MPIEXEC_NUMPROC_FLAG "-n" CACHE
            STRING "flag used to specify number of processes." FORCE)
      # else()
      endif()

      # Mark some of the variables created by the above logic as
      # 'advanced' so that they do not show up in the 'simple' ccmake 
      # view. 
      mark_as_advanced( MPI_EXTRA_LIBRARY MPI_LIBRARY )
      set( file_cmd ${file_cmd} CACHE INTERNAL "file command" )

      message(STATUS "Looking for MPI...found")

      # Sanity Checks for DRACO_C4==MPI
      if( "${MPI_CORES_PER_CPU}x" STREQUAL "x" )
         message( FATAL_ERROR "setupMPILibrariesUnix:: MPI_CORES_PER_CPU is not set!")
      endif()

   endif( NOT "${DRACO_C4}" STREQUAL "SCALAR" )

   set( MPI_SETUP_DONE ON CACHE INTERNAL "Have we completed the MPI setup call?" )

   # We use variables like ${MPI_CXX_INCLUDE_PATH} with the
   # assumption that the stirng will be empty if MPI is not found.
   # Make sure the string is empty...

   if( "${MPI_CXX_INCLUDE_PATH}" STREQUAL "MPI_CXX_INCLUDE_PATH-NOTFOUND")
      set( MPI_HEADER_PATH "" CACHE PATH 
         "MPI not found, empty value." FORCE )
      set( MPI_CXX_INCLUDE_PATH "" CACHE PATH 
         "MPI not found, empty value." FORCE )
      set( MPI_C_INCLUDE_PATH "" CACHE PATH 
         "MPI not found, empty value." FORCE )
   endif()

endmacro()

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
   foreach( version 3.4.1 3.4.2 )   
      if( EXISTS  ${LAPACK_LIB_DIR}/cmake/lapack-${version} )
         list( APPEND CMAKE_PREFIX_PATH ${LAPACK_LIB_DIR}/cmake/lapack-${version} )
         find_package( lapack CONFIG )
      endif()
   endforeach()
   foreach( config NOCONFIG DEBUG RELEASE )
      get_target_property(tmp lapack IMPORTED_LOCATION_${config} )
      if( EXISTS ${tmp} )
         set( lapack_FOUND TRUE )
      endif()
   endforeach()
   if( lapack_FOUND)
      message( STATUS "Looking for lapack...found")
      set( lapack_FOUND ${lapack_FOUND} CACHE BOOL "Did we find LAPACK." FORCE )
   else()
      message( "Looking for lapack...not found")
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
   mark_as_advanced( 
      CUDA_SDK_ROOT_DIR 
      CUDA_VERBOSE_BUILD
      CUDA_TOOLKIT_ROOT_DIR 
      CUDA_BUILD_CUBIN
      CUDA_BUILD_EMULATION
      )

endmacro()

#------------------------------------------------------------------------------
# Helper macros for setup_global_libraries()
#------------------------------------------------------------------------------
macro( SetupVendorLibrariesUnix )

   # GSL ----------------------------------------------------------------------
   # message( STATUS "Looking for GSL...")
   if( DRACO_LIBRARY_TYPE MATCHES "STATIC" )
      set( GSL_STATIC ON )
   endif()
   find_package( GSL QUIET )

   # Random123 ----------------------------------------------------------------
   # message( STATUS "Looking for Random123...")
   find_package( Random123 QUIET )

   # GRACE ------------------------------------------------------------------
   find_package( Grace QUIET )
   set_package_properties( Grace PROPERTIES
      DESCRIPTION "A WYSIWYG 2D plotting tool."
      TYPE OPTIONAL
      PURPOSE "Required for bulding the plot2D component."
      )

   # CUDA ------------------------------------------------------------------
   setupCudaEnv()
   
endmacro()

#------------------------------------------------------------------------------
# Helper macros for setup_global_libraries()
#------------------------------------------------------------------------------
macro( setupMPILibrariesWindows )

   # MPI ---------------------------------------------------------------------
   if( NOT "${DRACO_C4}" STREQUAL "SCALAR" ) # AND "${MPIEXEC}x" STREQUAL "x" )
      # Try to find MPI in the default locations (look for mpic++ in PATH)
      # This module will set the following variables:
      #   MPI_FOUND                  TRUE have we found MPI
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
      
      message(STATUS "Looking for MPI...")
      find_package( MPI )

      # Second chance using $MPIRUN (old Draco setup format -- ask JDD).
      if( NOT ${MPI_FOUND} AND EXISTS "${MPIRUN}" )
         set( MPIEXEC $ENV{MPIRUN} )
         find_package( MPI )
      endif()

      # Set Draco build system variables based on what we know about MPI.
      if( MPI_FOUND )
         set( DRACO_C4 "MPI" )  
         if( NOT MPIEXEC )
            message( FATAL_ERROR 
               "MPI found but mpirun not in PATH. Aborting" )
         endif()
         # Try to find the fortran mpi library
         if( EXISTS ${MPI_LIB_DIR} )
            find_library( MPI_Fortran_LIB mpi_f77 HINTS ${MPI_LIB_DIR} )
            mark_as_advanced( MPI_Fortran_LIB )
         endif()
         unset( C4_SCALAR )
      else()
         set( DRACO_C4 "SCALAR" )
         unset( C4_MPI )
         set( MPI_INCLUDE_PATH "" CACHE FILEPATH 
            "no mpi library for scalar build." FORCE )
         set( MPI_LIBRARY "" CACHE FILEPATH 
            "no mpi library for scalar build." FORCE )
         set( MPI_LIBRARIES "" CACHE FILEPATH 
            "no mpi library for scalar build." FORCE )
      endif()

      # Save the result in the cache file.
      set( DRACO_C4 "${DRACO_C4}" CACHE STRING 
         "C4 communication mode (SCALAR or MPI)" )
      if( "${DRACO_C4}" STREQUAL "MPI"    OR 
            "${DRACO_C4}" STREQUAL "SCALAR" )
      else()
         message( FATAL_ERROR "DRACO_C4 must be either MPI or SCALAR" )
      endif()
      
      # Check flavor and add optional flags
      if( "${MPIEXEC}" MATCHES openmpi )
         set( MPI_FLAVOR "openmpi" CACHE STRING "Flavor of MPI." )
         
         set( MPI_CORES_PER_CPU 4 )
         set( MPI_PHYSICAL_CORES 0 )
         math( EXPR MPI_CPUS_PER_NODE "${MPIEXEC_MAX_NUMPROCS} / ${MPI_CORES_PER_CPU}" )
         set( MPI_CPUS_PER_NODE ${MPI_CPUS_PER_NODE} CACHE STRING
            "Number of multi-core CPUs per node" FORCE )
         set( MPI_CORES_PER_CPU ${MPI_CORES_PER_CPU} CACHE STRING
            "Number of cores per cpu" FORCE )
            
         # Check for hyperthreading - This is important for reserving
         # threads for OpenMP tests...

         # correct base-zero indexing
         math( EXPR MPI_PHYSICAL_CORES "${MPI_PHYSICAL_CORES} + 1" )
         math( EXPR MPI_MAX_NUMPROCS_PHYSICAL
            "${MPI_PHYSICAL_CORES} * ${MPI_CORES_PER_CPU}" )
         if( "${MPI_MAX_NUMPROCS_PHYSICAL}" STREQUAL "${MPIEXEC_MAX_NUMPROCS}" )
            set( MPI_HYPERTHREADING "OFF" CACHE BOOL "Are we using hyperthreading?" FORCE )
         else()
            set( MPI_HYPERTHREADING "ON" CACHE BOOL "Are we using hyperthreading?" FORCE )
         endif()  
         
         #
         # EAP's flags can be found in Test.rh/General/run_job.pl
         # (look for $other_args).  In particular, it may be useful to
         # examine EAP's options for srun or aprun.
         # 
         # \sa
         # http://blogs.cisco.com/performance/open-mpi-v1-5-processor-affinity-options/
         #
         # --cpus-per-proc <N> option is needed for multi-threaded
         #   processes.  Use this as "threads per MPI rank."
         # --bind-to-socket will bind MPI ranks ordered by socket first (rank 0 onto
         #   socket 0, then rank 1 onto socket 1, etc.)

         # --bind-to-core added in OpenMPI-1.4
         set( MPIEXEC_OMP_POSTFLAGS 
            -bind-to-socket -cpus-per-proc ${MPI_CORES_PER_CPU} --report-bindings
            CACHE STRING "extra mpirun flags (list)." FORCE )
         set( MPIEXEC_OMP_POSTFLAGS_STRING 
            "-bind-to-socket -cpus-per-proc ${MPI_CORES_PER_CPU} --report-bindings"
            CACHE STRING "extra mpirun flags (list)." FORCE)
         mark_as_advanced( MPI_FLAVOR MPIEXEC_OMP_POSTFLAGS_STRING MPIEXEC_OMP_POSTFLAGS
            MPI_LIBRARIES )    
            
      elseif("${MPIEXEC}" MATCHES "Microsoft HPC" )
         set( MPI_FLAVOR "MicrosoftHPC" CACHE STRING "Flavor of MPI." )
         
         set( MPI_CORES_PER_CPU 4 )
         set( MPI_PHYSICAL_CORES 0 )
         math( EXPR MPI_CPUS_PER_NODE "${MPIEXEC_MAX_NUMPROCS} / ${MPI_CORES_PER_CPU}" )
         set( MPI_CPUS_PER_NODE ${MPI_CPUS_PER_NODE} CACHE STRING
            "Number of multi-core CPUs per node" FORCE )
         set( MPI_CORES_PER_CPU ${MPI_CORES_PER_CPU} CACHE STRING
            "Number of cores per cpu" FORCE )
            
         # Check for hyperthreading - This is important for reserving
         # threads for OpenMP tests...

         # correct base-zero indexing
         math( EXPR MPI_PHYSICAL_CORES "${MPI_PHYSICAL_CORES} + 1" )
         math( EXPR MPI_MAX_NUMPROCS_PHYSICAL
            "${MPI_PHYSICAL_CORES} * ${MPI_CORES_PER_CPU}" )
         if( "${MPI_MAX_NUMPROCS_PHYSICAL}" STREQUAL "${MPIEXEC_MAX_NUMPROCS}" )
            set( MPI_HYPERTHREADING "OFF" CACHE BOOL "Are we using hyperthreading?" FORCE )
         else()
            set( MPI_HYPERTHREADING "ON" CACHE BOOL "Are we using hyperthreading?" FORCE )
         endif()  
         
         #
         # EAP's flags can be found in Test.rh/General/run_job.pl
         # (look for $other_args).  In particular, it may be useful to
         # examine EAP's options for srun or aprun.
         # 
         # \sa
         # http://blogs.cisco.com/performance/open-mpi-v1-5-processor-affinity-options/
         #
         # --cpus-per-proc <N> option is needed for multi-threaded
         #   processes.  Use this as "threads per MPI rank."
         # --bind-to-socket will bind MPI ranks ordered by socket first (rank 0 onto
         #   socket 0, then rank 1 onto socket 1, etc.)

         # --bind-to-core added in OpenMPI-1.4
         set( MPIEXEC_OMP_POSTFLAGS 
            -bind-to-socket -cpus-per-proc ${MPI_CORES_PER_CPU} --report-bindings
            CACHE STRING "extra mpirun flags (list)." FORCE )
         set( MPIEXEC_OMP_POSTFLAGS_STRING 
            "-bind-to-socket -cpus-per-proc ${MPI_CORES_PER_CPU} --report-bindings"
            CACHE STRING "extra mpirun flags (list)." FORCE)
         mark_as_advanced( MPI_FLAVOR MPIEXEC_OMP_POSTFLAGS_STRING MPIEXEC_OMP_POSTFLAGS
            MPI_LIBRARIES )    
      endif()
      
   endif() # NOT "${DRACO_C4}" STREQUAL "SCALAR"

   set( MPI_SETUP_DONE ON CACHE INTERNAL "Have we completed the MPI setup call?" )

endmacro( setupMPILibrariesWindows )

macro( SetupVendorLibrariesWindows )

   # GSL ---------------------------------------------------------------------
   message( STATUS "Looking for GSL...")
   set( GSL_INC_DIR "${VENDOR_DIR}/gsl/include" )
   set( GSL_LIB_DIR "${VENDOR_DIR}/gsl/lib" )
   
   # Use static BLAS libraries
   set(GSL_STATIC ON)   
   find_package( GSL REQUIRED )
   
   # if( GSL_FOUND )
      # message( STATUS "Looking for GSL...   FOUND")
   # else()
      # message( STATUS "Looking for GSL...   NOT FOUND")
   # endif()

   # Random123 ---------------------------------------------------------------
   message( STATUS "Looking for Random123...")
   find_package( Random123 REQUIRED )

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

    if( NOT GSL_LIB_DIR AND IS_DIRECTORY $ENV{GSL_LIB_DIR}  )
        set( GSL_LIB_DIR $ENV{GSL_LIB_DIR} )
        set( GSL_INC_DIR $ENV{GSL_INC_DIR} )
    endif()
    if( NOT GSL_LIB_DIR AND IS_DIRECTORY ${VENDOR_DIR}/gsl/lib )
        set( GSL_LIB_DIR "${VENDOR_DIR}/gsl/lib" )
        set( GSL_INC_DIR "${VENDOR_DIR}/gsl/include" )
    endif()

    if( NOT RANDOM123_INC_DIR AND IS_DIRECTORY $ENV{RANDOM123_INC_DIR}  )
        set( RANDOM123_INC_DIR $ENV{RANDOM123_INC_DIR} )
    endif()
    if( NOT RANDOM123_INC_DIR AND IS_DIRECTORY ${VENDOR_DIR}/random123-1.06/include )
        set( RANDOM123_INC_DIR "${VENDOR_DIR}/random123-1.06/include" )
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
    set_package_properties( GSL PROPERTIES
       URL "http://www.gnu.org/s/gsl/"
       DESCRIPTION "GNU Scientific Library"
       TYPE REQUIRED
       PURPOSE "Required for bulding quadrature and rng components."
       )  
    set_package_properties( Random123 PROPERTIES
       URL "http://www.deshawresearch.com/resources_random123.html"
       DESCRIPTION "a library of counter-based random number generators"
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
