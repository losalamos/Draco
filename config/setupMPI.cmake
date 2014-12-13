#-----------------------------*-cmake-*----------------------------------------#
# file   config/setupMPI.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   2014 Sep 22
# brief  Setup MPI Vendors
# note   Copyright (C) 2014 Los Alamos National Security, LLC.
#        All rights reserved.
#
# Try to find MPI in the default locations (look for mpic++ in PATH)
# This module will set the following variables:
#
# MPI_<lang>_FOUND           TRUE if FindMPI found MPI flags for <lang>
# MPI_<lang>_COMPILER        MPI Compiler wrapper for <lang>
# MPI_<lang>_COMPILE_FLAGS   Compilation flags for MPI programs
# MPI_<lang>_INCLUDE_PATH    Include path(s) for MPI header
# MPI_<lang>_LINK_FLAGS      Linking flags for MPI programs
# MPI_<lang>_LIBRARIES       All libraries to link MPI programs against
#
# MPIEXEC                    Executable for running MPI programs
# MPIEXEC_NUMPROC_FLAG       Flag to pass to MPIEXEC before giving it the
#                            number of processors to run on
# MPIEXEC_PREFLAGS           Flags to pass to MPIEXEC directly before the
#                            executable to run.
# MPIEXEC_POSTFLAGS          Flags to pass to MPIEXEC after all other flags.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

include( FeatureSummary )

##---------------------------------------------------------------------------##
## Set Draco specific MPI variables
##---------------------------------------------------------------------------##
macro( setupDracoMPIVars )

      # Set Draco build system variables based on what we know about MPI.
      if( MPI_FOUND )
         set( DRACO_C4 "MPI" )
      else()
         set( DRACO_C4 "SCALAR" )
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

      # Find the version
      execute_process( COMMAND ${MPIEXEC} --version
         OUTPUT_VARIABLE DBS_MPI_VER_OUT
         ERROR_VARIABLE DBS_MPI_VER_ERR)

      set( DBS_MPI_VER "${DBS_MPI_VER_OUT}${DBS_MPI_VER_ERR}")

endmacro()

##---------------------------------------------------------------------------##
## Query CPU topology
##
## Returns:
##   MPI_CORES_PER_CPU
##   MPI_CPUS_PER_NODE
##   MPI_PHYSICAL_CORES
##   MPI_MAX_NUMPROCS_PHYSICAL
##   MPI_HYPERTHREADING
##
## See also:
##   - Try running 'lstopo' for a graphical view of the local
##     topology.
##   - EAP's flags can be found in Test.rh/General/run_job.pl (look
##     for $other_args).  In particular, it may be useful to examine
##     EAP's options for srun or aprun.
##   - http://blogs.cisco.com/performance/open-mpi-v1-5-processor-affinity-options/
##
##---------------------------------------------------------------------------##
macro( query_topology )

  # start with default values
  set( MPI_CORES_PER_CPU 4 )
  set( MPI_PHYSICAL_CORES 0 )

  # read the system's cpuinfo...
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

  math( EXPR MPI_CPUS_PER_NODE
    "${MPIEXEC_MAX_NUMPROCS} / ${MPI_CORES_PER_CPU}" )
  set( MPI_CPUS_PER_NODE ${MPI_CPUS_PER_NODE} CACHE STRING
    "Number of multi-core CPUs per node" FORCE )
  set( MPI_CORES_PER_CPU ${MPI_CORES_PER_CPU} CACHE STRING
    "Number of cores per cpu" FORCE )

  #
  # Check for hyperthreading - This is important for reserving threads for OpenMP tests...
  #
  # correct base-zero indexing
  math( EXPR MPI_PHYSICAL_CORES        "${MPI_PHYSICAL_CORES} + 1" )
  math( EXPR MPI_MAX_NUMPROCS_PHYSICAL "${MPI_PHYSICAL_CORES} * ${MPI_CORES_PER_CPU}" )
  if( "${MPI_MAX_NUMPROCS_PHYSICAL}" STREQUAL "${MPIEXEC_MAX_NUMPROCS}" )
    set( MPI_HYPERTHREADING "OFF" CACHE BOOL "Are we using hyperthreading?" FORCE )
  else()
    set( MPI_HYPERTHREADING "ON" CACHE BOOL "Are we using hyperthreading?" FORCE )
  endif()

  #
  # Setup for OMP plus MPI
  #
  if( MPI_FLAVOR STREQUAL "openmpi" )
    if( ${DBS_MPI_VER_MAJOR}.${DBS_MPI_VER_MINOR} VERSION_LESS 1.7 )
      set( MPIEXEC_OMP_POSTFLAGS
        "-bind-to-socket -cpus-per-proc ${MPI_CORES_PER_CPU} --report-bindings" )
    elseif( ${DBS_MPI_VER_MAJOR}.${DBS_MPI_VER_MINOR} VERSION_GREATER 1.7 )
      set( MPIEXEC_OMP_POSTFLAGS
        "-bind-to socket --map-by ppr:${MPI_CORES_PER_CPU}:socket --report-bindings" )
    else()  # Version 1.7.4
      set( MPIEXEC_OMP_POSTFLAGS
        "-bind-to socket --map-by socket:PPR=${MPI_CORES_PER_CPU}" )
    endif()

    set( MPIEXEC_OMP_POSTFLAGS ${MPIEXEC_OMP_POSTFLAGS}
      CACHE STRING "extra mpirun flags (list)." FORCE )
  endif()

  mark_as_advanced( MPI_CPUS_PER_NODE MPI_CORES_PER_CPU
     MPI_PHYSICAL_CORES MPI_MAX_NUMPROCS_PHYSICAL MPI_HYPERTHREADING )

endmacro()

#------------------------------------------------------------------------------#
# Setup MPI when on Linux
#------------------------------------------------------------------------------#
macro( setupMPILibrariesUnix )

   # MPI ---------------------------------------------------------------------
   if( NOT "${DRACO_C4}" STREQUAL "SCALAR" )

      message(STATUS "Looking for MPI...")

      # Preserve data that may already be set.
      if( DEFINED ENV{MPIEXEC} )
        set( MPIEXEC $ENV{MPIEXEC} )
      endif()

      # Temporary work around until FindMPI.cmake is fixed:
      # Setting MPI_<LANG>_COMPILER and MPI_<LANG>_NO_INTERROGATE
      # forces FindMPI to skip it's bad logic and just rely on the MPI
      # compiler wrapper to do the right thing. see Bug #467.
      foreach( lang C CXX Fortran )
        if( "${CMAKE_${lang}_COMPILER}" MATCHES "mpi[A-z+]+" )
          get_filename_component( compiler_wo_path "${CMAKE_${lang}_COMPILER}" NAME )
          set( MPI_${lang}_COMPILER ${CMAKE_${lang}_COMPILER} )
          set( MPI_${lang}_NO_INTERROGATE ${CMAKE_${lang}_COMPILER} )
        endif()
      endforeach()

      # Call the standard CMake FindMPI macro.
      find_package( MPI QUIET )

      # Set DRACO_C4 and other variables
      setupDracoMPIVars()

      # -------------------------------------------------------------------------------- #
      # Check flavor and add optional flags
      if( "${MPIEXEC}"     MATCHES openmpi   OR
          "${DBS_MPI_VER}" MATCHES open-mpi )

        set( MPI_FLAVOR "openmpi" CACHE STRING "Flavor of MPI." ) # OpenMPI

         # Find the version of OpenMPI

         if( "${DBS_MPI_VER}" MATCHES "[0-9][.][0-9][.][0-9]" )
            string( REGEX REPLACE ".*([0-9])[.]([0-9])[.]([0-9]).*" "\\1"
               DBS_MPI_VER_MAJOR ${DBS_MPI_VER} )
            string( REGEX REPLACE ".*([0-9])[.]([0-9])[.]([0-9]).*" "\\2"
               DBS_MPI_VER_MINOR ${DBS_MPI_VER} )
         elseif( "${DBS_MPI_VER}" MATCHES "[0-9][.][0-9]" )
            string( REGEX REPLACE ".*([0-9])[.]([0-9]).*" "\\1"
               DBS_MPI_VER_MAJOR ${DBS_MPI_VER} )
            string( REGEX REPLACE ".*([0-9])[.]([0-9]).*" "\\2"
               DBS_MPI_VER_MINOR ${DBS_MPI_VER} )
         endif()

         # sanity check, these OpenMPI flags (below) require version >= 1.4
         if( ${DBS_MPI_VER_MAJOR}.${DBS_MPI_VER_MINOR} VERSION_LESS 1.4 )
            message( FATAL_ERROR "OpenMPI version < 1.4 found." )
         endif()

         # Setting mpi_paffinity_alone to 0 allows parallel ctest to
         # work correctly.  MPIEXEC_POSTFLAGS only affects MPI-only
         # tests (and not MPI+OpenMP tests).

         # This flag also shows up in
         # jayenne/pkg_tools/run_milagro_test.py and regress_funcs.py.
         if( ${DBS_MPI_VER_MAJOR}.${DBS_MPI_VER_MINOR} VERSION_LESS 1.7 )
           set( MPIEXEC_POSTFLAGS "--mca mpi_paffinity_alone 0" CACHE
               STRING "extra mpirun flags (list)." FORCE)
         else()
           # Flag provided by Sam Gutierrez (2014-04-08),
           set( MPIEXEC_POSTFLAGS "-mca hwloc_base_binding_policy none" CACHE
               STRING "extra mpirun flags (list)." FORCE)
         endif()

         # Find cores/cpu, cpu/node, hyperthreading
         query_topology()

         #
         # Setup for OMP plus MPI
         #

         if( ${DBS_MPI_VER_MAJOR}.${DBS_MPI_VER_MINOR} VERSION_LESS 1.7 )
           set( MPIEXEC_OMP_POSTFLAGS
             "-bind-to-socket -cpus-per-proc ${MPI_CORES_PER_CPU} --report-bindings"
             CACHE STRING "extra mpirun flags (list)." FORCE )
         elseif( ${DBS_MPI_VER_MAJOR}.${DBS_MPI_VER_MINOR} VERSION_GREATER 1.7 )
           set( MPIEXEC_OMP_POSTFLAGS
             "-bind-to socket --map-by ppr:${MPI_CORES_PER_CPU}:socket --report-bindings"
             CACHE STRING "extra mpirun flags (list)." FORCE )
        else()
          # Version 1.7.4
           set( MPIEXEC_OMP_POSTFLAGS
             "-bind-to socket --map-by socket:PPR=${MPI_CORES_PER_CPU}"
             CACHE STRING "extra mpirun flags (list)." FORCE )
         endif()

         mark_as_advanced( MPI_CORES_PER_CPU MPI_PHYSICAL_CORES
           MPI_CPUS_PER_NODE MPI_MAX_NUMPROCS_PHYSICAL MPI_HYPERTHREADING
           MPIEXEC_OMP_POSTFLAGS )

         mark_as_advanced( MPI_FLAVOR MPIEXEC_OMP_POSTFLAGS MPI_LIBRARIES )

      ### Cray wrappers for mpirun
      elseif( "${MPIEXEC}" MATCHES aprun)
         # According to email from Mike McKay (2013/04/18), we might
         # need to set the the mpirun command to something like:
         #
         #   setenv OMP_NUM_THREADS ${MPI_CORES_PER_CPU}; aprun ... -d ${MPI_CORES_PER_CPU}
         #
         # consider '-cc none'
         if( "${DBS_MPI_VER}x" STREQUAL "x" AND
             NOT "$ENV{CRAY_MPICH2_VER}x" STREQUAL "x" )
           set( DBS_MPI_VER $ENV{CRAY_MPICH2_VER} )
           if( "${DBS_MPI_VER}" MATCHES "[0-9][.][0-9][.][0-9]" )
            string( REGEX REPLACE ".*([0-9])[.]([0-9])[.]([0-9]).*" "\\1"
               DBS_MPI_VER_MAJOR ${DBS_MPI_VER} )
            string( REGEX REPLACE ".*([0-9])[.]([0-9])[.]([0-9]).*" "\\2"
               DBS_MPI_VER_MINOR ${DBS_MPI_VER} )
             endif()
         endif()
         if( NOT "$ENV{OMP_NUM_THREADS}x" STREQUAL "x" )
            set( MPI_CPUS_PER_NODE 1 CACHE STRING
               "Number of multi-core CPUs per node" FORCE )
            set( MPI_CORES_PER_CPU $ENV{OMP_NUM_THREADS} CACHE STRING
               "Number of cores per cpu" FORCE )
            set( MPIEXEC_OMP_POSTFLAGS "-d ${MPI_CORES_PER_CPU}" CACHE
               STRING "extra mpirun flags (list)." FORCE)
         else()
            message( STATUS "
WARNING: ENV{OMP_NUM_THREADS} is not set in your environment,
         all OMP tests will be disabled." )
         endif()

      ### Sequoia
      elseif( "${MPIEXEC}" MATCHES srun)
         if( NOT "$ENV{OMP_NUM_THREADS}x" STREQUAL "x" )
            set( MPI_CPUS_PER_NODE 1 CACHE STRING
               "Number of multi-core CPUs per node" FORCE )
            set( MPI_CORES_PER_CPU $ENV{OMP_NUM_THREADS} CACHE STRING
               "Number of cores per cpu" FORCE )
            set( MPIEXEC_OMP_POSTFLAGS "-c${MPI_CORES_PER_CPU}" CACHE
               STRING "extra mpirun flags (list)." FORCE)
         else()
            message( STATUS "
WARNING: ENV{OMP_NUM_THREADS} is not set in your environment,
         all OMP tests will be disabled." )
         endif()

         set( MPIEXEC_NUMPROC_FLAG "-n" CACHE
           STRING "flag used to specify number of processes." FORCE)

      ### Intel MPI
      elseif( "${MPIEXEC}" MATCHES intel-mpi OR
          "${DBS_MPI_VER}" MATCHES "Intel[(]R[)] MPI Library" )

         set( MPI_FLAVOR "impi" CACHE STRING "Flavor of MPI." )

         # Find the version of Intel MPI

         if( "${DBS_MPI_VER}" MATCHES "[0-9][.][0-9][.][0-9]" )
            string( REGEX REPLACE ".*([0-9])[.]([0-9])[.]([0-9]).*" "\\1"
               DBS_MPI_VER_MAJOR ${DBS_MPI_VER} )
            string( REGEX REPLACE ".*([0-9])[.]([0-9])[.]([0-9]).*" "\\2"
               DBS_MPI_VER_MINOR ${DBS_MPI_VER} )
         elseif( "${DBS_MPI_VER}" MATCHES "[0-9][.][0-9]" )
            string( REGEX REPLACE ".*([0-9])[.]([0-9]).*" "\\1"
               DBS_MPI_VER_MAJOR ${DBS_MPI_VER} )
            string( REGEX REPLACE ".*([0-9])[.]([0-9]).*" "\\2"
               DBS_MPI_VER_MINOR ${DBS_MPI_VER} )
         endif()

         # Find cores/cpu, cpu/node, hyperthreading
         query_topology()

      else()
         message( FATAL_ERROR "
The Draco build system doesn't know how to configure the build for
MPIEXEC=${MPIEXEC}")
      endif()

      # Mark some of the variables created by the above logic as
      # 'advanced' so that they do not show up in the 'simple' ccmake
      # view.
      mark_as_advanced( MPI_EXTRA_LIBRARY MPI_LIBRARY )
      set( file_cmd ${file_cmd} CACHE INTERNAL "file command" )

      message(STATUS "Looking for MPI.......found ${MPIEXEC}")

      # Sanity Checks for DRACO_C4==MPI
      if( "${MPI_CORES_PER_CPU}x" STREQUAL "x" )
         message( FATAL_ERROR "setupMPILibrariesUnix:: MPI_CORES_PER_CPU is not set!")
      endif()

   endif( NOT "${DRACO_C4}" STREQUAL "SCALAR" )

   set( MPI_SETUP_DONE ON CACHE INTERNAL "Have we completed the MPI setup call?" )


endmacro()

##---------------------------------------------------------------------------##
##
##---------------------------------------------------------------------------##

macro( setupMPILibrariesWindows )

   # MPI ---------------------------------------------------------------------
   if( NOT "${DRACO_C4}" STREQUAL "SCALAR" )

      message(STATUS "Looking for MPI...")
      find_package( MPI )

      # Second chance using $MPIRUN (old Draco setup format -- ask JDD).
      if( NOT ${MPI_FOUND} AND EXISTS "${MPIRUN}" )
         set( MPIEXEC $ENV{MPIRUN} )
         find_package( MPI )
      endif()

      setupDracoMPIVars()

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
         if( ${DBS_MPI_VER_MAJOR}.${DBS_MPI_VER_MINOR} VERSION_LESS 1.7 )
           set( MPIEXEC_OMP_POSTFLAGS
             "-bind-to-socket -cpus-per-proc ${MPI_CORES_PER_CPU} --report-bindings"
             CACHE STRING "extra mpirun flags (list)." FORCE )
         else()
           set( MPIEXEC_OMP_POSTFLAGS
             "-bind-to socket --map-by ppr:${MPI_CORES_PER_CPU}:socket --report-bindings"
             CACHE STRING "extra mpirun flags (list)." FORCE )
         endif()
         mark_as_advanced( MPI_FLAVOR MPIEXEC_OMP_POSTFLAGS MPI_LIBRARIES )

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
         set( MPIEXEC_OMP_POSTFLAGS "-exitcodes"
            CACHE STRING "extra mpirun flags (list)." FORCE )
         mark_as_advanced( MPI_FLAVOR MPIEXEC_OMP_POSTFLAGS MPI_LIBRARIES )
      elseif("${MPIEXEC}" MATCHES "MPICH2" )
         set( MPI_FLAVOR "MPICH2" CACHE STRING "Flavor of MPI." )

         include(ProcessorCount)
         ProcessorCount(MPIEXEC_MAX_NUMPROCS)
         set( MPI_CORES_PER_CPU ${MPIEXEC_MAX_NUMPROCS} )
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
         set( MPIEXEC_OMP_POSTFLAGS "-exitcodes"
            CACHE STRING "extra mpirun flags (list)." FORCE )
         mark_as_advanced( MPI_FLAVOR MPIEXEC_OMP_POSTFLAGS MPI_LIBRARIES )
      endif()

   endif() # NOT "${DRACO_C4}" STREQUAL "SCALAR"

   set( MPI_SETUP_DONE ON CACHE INTERNAL "Have we completed the MPI setup call?" )
   if( ${MPI_FOUND} )
      message(STATUS "Looking for MPI...${MPIEXEC}")
   else()
      message(STATUS "Looking for MPI...not found")
   endif()

endmacro( setupMPILibrariesWindows )

#----------------------------------------------------------------------#
# End setupMPI.cmake
#----------------------------------------------------------------------#
