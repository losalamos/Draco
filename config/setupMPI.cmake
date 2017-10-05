#-----------------------------*-cmake-*----------------------------------------#
# file   config/setupMPI.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   2016 Sep 22
# brief  Setup MPI Vendors
# note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
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
#
# DRACO_C4                   MPI|SCALAR
# C4_SCALAR                  BOOL
# C4_MPI                     BOOL
#
#------------------------------------------------------------------------------#
include( FeatureSummary )

##---------------------------------------------------------------------------##
## Set Draco specific MPI variables
##---------------------------------------------------------------------------##
macro( setupDracoMPIVars )

  # Set Draco build system variables based on what we know about MPI.
  if( MPI_CXX_FOUND )
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

  if( "${DRACO_C4}" MATCHES "MPI" )
    set( C4_SCALAR 0 )
    set( C4_MPI 1 )
  else()
    set( C4_SCALAR 1 )
    set( C4_MPI 0 )
  endif()
  set( C4_SCALAR ${C4_SCALAR} CACHE STRING
    "Are we building a scalar-only version (no mpi in c4/config.h)?"
    FORCE )
  set( C4_MPI ${C4_MPI} CACHE STRING
    "Are we building an MPI aware version? (c4/config.h)" FORCE )
  mark_as_advanced( C4_MPI C4_SCALAR )

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
##   - Try running 'lstopo' for a graphical view of the local topology or
##     'lscpu' for a text version.
##   - EAP's flags can be found in Test.rh/General/run_job.pl (look for
##     $other_args).  In particular, it may be useful to examine EAP's options
##     for srun or aprun.
##---------------------------------------------------------------------------##
macro( query_topology )

  # start with default values
  set( MPI_CORES_PER_CPU 4 )
  set( MPI_PHYSICAL_CORES 1 )

  if( "${SITENAME}" STREQUAL "Trinitite" OR
      "${SITENAME}" STREQUAL "Trinity" )
    # Backend is different than build-node
    if( $ENV{CRAY_CPU_TARGET} MATCHES "mic-knl" )
      set( MPI_CORES_PER_CPU 17 )
      set( MPI_PHYSICAL_CORES 4 )
      set( MPIEXEC_MAX_NUMPROCS 68 CACHE STRING "Max procs on node." FORCE )
    else()
      set( MPI_CORES_PER_CPU 16 )
      set( MPI_PHYSICAL_CORES 2 )
      set( MPIEXEC_MAX_NUMPROCS 32 CACHE STRING "Max procs on node." FORCE )
    endif()
  elseif( EXISTS "/proc/cpuinfo" )
    # read the system's cpuinfo...
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
    # correct 0-based indexing
    math( EXPR MPI_PHYSICAL_CORES "${MPI_PHYSICAL_CORES} + 1" )
  endif()

  math( EXPR MPI_CPUS_PER_NODE
    "${MPIEXEC_MAX_NUMPROCS} / ${MPI_CORES_PER_CPU}" )
  set( MPI_CPUS_PER_NODE ${MPI_CPUS_PER_NODE} CACHE STRING
    "Number of multi-core CPUs per node" FORCE )
  set( MPI_CORES_PER_CPU ${MPI_CORES_PER_CPU} CACHE STRING
    "Number of cores per cpu" FORCE )

  #
  # Check for hyperthreading - This is important for reserving threads for
  # OpenMP tests...
  #
  math( EXPR MPI_MAX_NUMPROCS_PHYSICAL "${MPI_PHYSICAL_CORES} * ${MPI_CORES_PER_CPU}" )
  if( "${MPI_MAX_NUMPROCS_PHYSICAL}" STREQUAL "${MPIEXEC_MAX_NUMPROCS}" )
    set( MPI_HYPERTHREADING "OFF" CACHE BOOL "Are we using hyperthreading?" FORCE )
  else()
    set( MPI_HYPERTHREADING "ON" CACHE BOOL "Are we using hyperthreading?" FORCE )
  endif()

endmacro()

##---------------------------------------------------------------------------##
## Setup OpenMPI
##---------------------------------------------------------------------------##
macro( setupOpenMPI )

  set( MPI_FLAVOR "openmpi" CACHE STRING "Flavor of MPI." ) # OpenMPI

  # Find the version of OpenMPI

  if( "${DBS_MPI_VER}" MATCHES "[0-9]+[.][0-9]+[.][0-9]+" )
    string( REGEX REPLACE ".*([0-9]+)[.]([0-9]+)[.]([0-9]+).*" "\\1"
      DBS_MPI_VER_MAJOR ${DBS_MPI_VER} )
    string( REGEX REPLACE ".*([0-9]+)[.]([0-9]+)[.]([0-9]+).*" "\\2"
      DBS_MPI_VER_MINOR ${DBS_MPI_VER} )
  elseif( "${DBS_MPI_VER}" MATCHES "[0-9]+[.][0-9]+" )
    string( REGEX REPLACE ".*([0-9]+)[.]([0-9]+).*" "\\1"
      DBS_MPI_VER_MAJOR ${DBS_MPI_VER} )
    string( REGEX REPLACE ".*([0-9]+)[.]([0-9]+).*" "\\2"
      DBS_MPI_VER_MINOR ${DBS_MPI_VER} )
  endif()

  # sanity check, these OpenMPI flags (below) require version >= 1.4
  if( ${DBS_MPI_VER_MAJOR}.${DBS_MPI_VER_MINOR} VERSION_LESS 1.8 )
    message( FATAL_ERROR "OpenMPI version < 1.8 found." )
  endif()

  # Setting mpi_paffinity_alone to 0 allows parallel ctest to work correctly.
  # MPIEXEC_POSTFLAGS only affects MPI-only tests (and not MPI+OpenMP tests).
  if( "$ENV{GITLAB_CI}" STREQUAL "true" )
    set(runasroot "--allow-run-as-root")
  endif()

  # This flag also shows up in jayenne/pkg_tools/run_milagro_test.py and
  # regress_funcs.py.

  # (2017-01-13) Bugs in openmpi-1.10.x are mostly fixed. Remove flags used
  # to work around bugs: '-mca btl self,vader -mca timer_require_monotonic 0'
  set( MPIEXEC_POSTFLAGS "-bind-to none ${runasroot}" CACHE STRING
    "extra mpirun flags (list)." FORCE)

  # Find cores/cpu, cpu/node, hyperthreading
  query_topology()

  #
  # Setup for OMP plus MPI
  #
  if( NOT APPLE )
    # -bind-to fails on OSX, See #691
    set( MPIEXEC_OMP_POSTFLAGS
      "--map-by ppr:${MPI_CORES_PER_CPU}:socket --report-bindings ${runasroot}" )
      # "-bind-to socket --map-by ppr:${MPI_CORES_PER_CPU}:socket --report-bindings ${runasroot}"
  endif()

  set( MPIEXEC_OMP_POSTFLAGS ${MPIEXEC_OMP_POSTFLAGS}
    CACHE STRING "extra mpirun flags (list)." FORCE )

  mark_as_advanced( MPI_CPUS_PER_NODE MPI_CORES_PER_CPU
     MPI_PHYSICAL_CORES MPI_MAX_NUMPROCS_PHYSICAL MPI_HYPERTHREADING )

endmacro()

##---------------------------------------------------------------------------##
## Setup Intel MPI
##---------------------------------------------------------------------------##
macro( setupIntelMPI )

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

endmacro()

##---------------------------------------------------------------------------##
## Setup Cray MPI wrappers (APRUN)
##---------------------------------------------------------------------------##
macro( setupCrayMPI )

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

  query_topology()

  # salloc/sbatch options:
  # --------------------
  # -N        limit job to a single node.
  # --gres=craynetwork:0 This option allows more than one srun to be running at
  #           the same time on the Cray. There are 4 gres “tokens” available. If
  #           unspecified, each srun invocation will consume all of
  #           them. Setting the value to 0 means consume none and allow the user
  #           to run as many concurrent jobs as there are cores available on the
  #           node. This should only be specified on the salloc/sbatch command.
  #           Gabe doesn't recommend this option for regression testing.
  # --vm-overcommit=disable|enable Do not allow overcommit of heap resources.
  # -p knl    Limit allocation to KNL nodes.
  # srun options:
  # --------------------
  # --cpu_bind=verbose,cores
  #           bind MPI ranks to cores
  #           print a summary of binding when run
  # --exclusive This option will keep concurrent jobs from running on the same
  #           cores. If you want to background tasks to have them run
  #           simultaneously, this option is required to be set or they will
  #           stomp on the same cores.

  set(postflags "-N 1 --cpu_bind=verbose,cores ")
  string(APPEND postflags " --exclusive")
  set( MPIEXEC_POSTFLAGS ${postflags} CACHE STRING
    "extra mpirun flags (list)." FORCE)

  set( MPIEXEC_OMP_POSTFLAGS "${MPIEXEC_POSTFLAGS} -c ${MPI_CORES_PER_CPU}"
    CACHE STRING "extra mpirun flags (list)." FORCE)

endmacro()

##---------------------------------------------------------------------------##
## Setup Sequoia MPI wrappers
##---------------------------------------------------------------------------##
macro( setupSequoiaMPI )

  if( NOT "$ENV{OMP_NUM_THREADS}x" STREQUAL "x" )
    set( MPI_CPUS_PER_NODE 1 CACHE STRING
      "Number of multi-core CPUs per node" FORCE )
    set( MPI_CORES_PER_CPU $ENV{OMP_NUM_THREADS} CACHE STRING
      "Number of cores per cpu" FORCE )
    set( MPIEXEC_OMP_POSTFLAGS "-c ${MPI_CORES_PER_CPU}" CACHE
      STRING "extra mpirun flags (list)." FORCE)
  else()
    message( STATUS "
WARNING: ENV{OMP_NUM_THREADS} is not set in your environment,
         all OMP tests will be disabled." )
  endif()

  set( MPIEXEC_NUMPROC_FLAG "-n" CACHE
    STRING "flag used to specify number of processes." FORCE)

  # [2016-11-17 KT] Ensure that the MPI headers in mpi.h use function signatures
  # from the MPI-2.2 standard that include 'const' parameters for MPI_put. The
  # MPI headers found on Sequoia claim to be MPI-2.2 compliant. However, this
  # CPP macro has been altered to be empty and this violates the v. 2.2
  # standard. For compatibility with the CCS-2 codebase, manually set this CPP
  # macro back to 'const'.  This definition will appear in c4/config.h.
  set( MPICH2_CONST "const" "Sequoia MPICH2-1.5 compile option.")
  mark_as_advanced( MPICH2_CONST )

endmacro()

#------------------------------------------------------------------------------#
# Setup MPI when on Linux
#------------------------------------------------------------------------------#
macro( setupMPILibrariesUnix )

   # MPI ---------------------------------------------------------------------
   if( NOT "${DRACO_C4}" STREQUAL "SCALAR" )

      message(STATUS "Looking for MPI...")

      # If this is a Cray system and the Cray MPI compile wrappers are used,
      # then do some special setup:

      if( CRAY_PE )
        set( MPIEXEC "srun" CACHE STRING
          "Program to execute MPI prallel programs." )
        set( MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING
          "mpirun flag used to specify the number of processors to use")
      endif()

      # Preserve data that may already be set.
      if( DEFINED ENV{MPIRUN} )
        set( MPIEXEC $ENV{MPIRUN} CACHE STRING
          "Program to execute MPI prallel programs." )
      elseif( DEFINED ENV{MPIEXEC} )
        set( MPIEXEC $ENV{MPIEXEC} CACHE STRING
          "Program to execute MPI prallel programs." )
      endif()

      # Temporary work around until FindMPI.cmake is fixed: Setting
      # MPI_<LANG>_COMPILER and MPI_<LANG>_NO_INTERROGATE forces FindMPI to skip
      # it's bad logic and just rely on the MPI compiler wrapper to do the right
      # thing. see Bug #467.
      foreach( lang C CXX Fortran )
        get_filename_component( CMAKE_${lang}_COMPILER_NOPATH
          "${CMAKE_${lang}_COMPILER}" NAME )
        if( "${CMAKE_${lang}_COMPILER_NOPATH}" MATCHES "^mpi[A-z+]+" )
          get_filename_component( compiler_wo_path "${CMAKE_${lang}_COMPILER}"
            NAME )
          set( MPI_${lang}_COMPILER ${CMAKE_${lang}_COMPILER} )
          set( MPI_${lang}_NO_INTERROGATE ${CMAKE_${lang}_COMPILER} )
        endif()
      endforeach()

      # Call the standard CMake FindMPI macro.
      find_package( MPI QUIET )

      # Set DRACO_C4 and other variables
      setupDracoMPIVars()

      # Find the mpirun version (skip this on Cray because this command seems to
      # occasionally hang).
      if( NOT CRAY_PE )
        execute_process( COMMAND ${MPIEXEC} --version
          OUTPUT_VARIABLE DBS_MPI_VER_OUT
          ERROR_VARIABLE DBS_MPI_VER_ERR)
        set( DBS_MPI_VER "${DBS_MPI_VER_OUT}${DBS_MPI_VER_ERR}")
      endif()

      set_package_properties( MPI PROPERTIES
        URL "http://www.open-mpi.org/"
        DESCRIPTION "A High Performance Message Passing Library"
        TYPE RECOMMENDED
        PURPOSE "If not available, all Draco components will be built as scalar applications."
        )

      # ---------------------------------------------------------------------- #
      # Check flavor and add optional flags
      #
      # Notes:
      # 1. For Intel MPI when cross compiling for MIC, the variable DBS_MPI_VER
      #    will be bad because this configuration is done on x86 but mpiexec is
      #    a mic executable and cannot be run on the x86.  To correctly match
      #    the MPI flavor to Intel (on darwin), we rely on the path MPIEXEC to
      #    match the string "impi/[0-9.]+/mic".

      if( "${MPIEXEC}" MATCHES openmpi OR "${DBS_MPI_VER}" MATCHES open-mpi )
        setupOpenMPI()

      elseif( "${MPIEXEC}" MATCHES intel-mpi OR
          "${MPIEXEC}" MATCHES "impi/.*/mic" OR
          "${DBS_MPI_VER}" MATCHES "Intel[(]R[)] MPI Library" )
        setupIntelMPI()

      elseif( CRAY_PE )
        setupCrayMPI()

      # LANL Cray systems also use srun, so this 'elseif' must appear after our
      # test for CRAY_PE.
      elseif( "${MPIEXEC}" MATCHES srun)
        setupSequoiaMPI()

      else()
         message( FATAL_ERROR "
The Draco build system doesn't know how to configure the build for
  MPIEXEC     = ${MPIEXEC}
  DBS_MPI_VER = ${DBS_MPI_VER}")
      endif()

      # Mark some of the variables created by the above logic as 'advanced' so
      # that they do not show up in the 'simple' ccmake view.
      mark_as_advanced( MPI_EXTRA_LIBRARY MPI_LIBRARY )

      message(STATUS "Looking for MPI.......found ${MPIEXEC}")

      # Sanity Checks for DRACO_C4==MPI
      if( "${MPI_CORES_PER_CPU}x" STREQUAL "x" )
         message( FATAL_ERROR "setupMPILibrariesUnix:: MPI_CORES_PER_CPU "
           "is not set!")
      endif()

    else()
      # Set DRACO_C4 and other variables
      setupDracoMPIVars()
    endif()

   mark_as_advanced( MPI_FLAVOR MPIEXEC_OMP_POSTFLAGS MPI_LIBRARIES )

 endmacro()

##---------------------------------------------------------------------------##
## setupMPILibrariesWindows
##---------------------------------------------------------------------------##

macro( setupMPILibrariesWindows )

   # MPI ---------------------------------------------------------------------
   if( NOT "${DRACO_C4}" STREQUAL "SCALAR" )

      message(STATUS "Looking for MPI...")
      find_package( MPI )

      # For MS-MPI 5, mpifptr.h is architecture dependent. Figure out
      # what arch this is and save this path to MPI_Fortran_INCLUDE_PATH
      list( GET MPI_CXX_LIBRARIES 0 first_cxx_mpi_library )
      if( first_cxx_mpi_library AND NOT MPI_Fortran_INCLUDE_PATH )
        get_filename_component( MPI_Fortran_INCLUDE_PATH "${first_cxx_mpi_library}" DIRECTORY )
        string( REPLACE "lib" "Include" MPI_Fortran_INCLUDE_PATH ${MPI_Fortran_INCLUDE_PATH} )
        set( MPI_Fortran_INCLUDE_PATH
             "${MPI_CXX_INCLUDE_PATH};${MPI_Fortran_INCLUDE_PATH}"
             CACHE STRING "Location for MPI include files for Fortran.")
      endif()

      setupDracoMPIVars()

      # Find the version
      # This is not working (hardwire it for now)
      execute_process( COMMAND "${MPIEXEC}" -help
        OUTPUT_VARIABLE DBS_MPI_VER_OUT
        ERROR_VARIABLE DBS_MPI_VER_ERR
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_STRIP_TRAILING_WHITESPACE
        )
      if( "${DBS_MPI_VER_OUT}" MATCHES "Microsoft MPI Startup Program" )
          string( REGEX REPLACE ".*Version ([0-9.]+).*" "\\1" DBS_MPI_VER "${DBS_MPI_VER_OUT}${DBS_MPI_VER_ERR}")
          string( REGEX REPLACE ".*([0-9])[.]([0-9])[.]([0-9]+).*" "\\1" DBS_MPI_VER_MAJOR ${DBS_MPI_VER} )
          string( REGEX REPLACE ".*([0-9])[.]([0-9])[.]([0-9]+).*" "\\2" DBS_MPI_VER_MINOR ${DBS_MPI_VER} )
          set( DBS_MPI_VER "${DBS_MPI_VER_MAJOR}.${DBS_MPI_VER_MINOR}")
      else()
         set(DBS_MPI_VER "5.0")
      endif()

      set_package_properties( MPI PROPERTIES
        URL "https://msdn.microsoft.com/en-us/library/bb524831%28v=vs.85%29.aspx"
        DESCRIPTION "Microsoft MPI"
        TYPE RECOMMENDED
        PURPOSE "If not available, all Draco components will be built as scalar applications."
        )

      # Check flavor and add optional flags
      if("${MPIEXEC}" MATCHES "Microsoft HPC" OR "${MPIEXEC}" MATCHES "Microsoft MPI")
         set( MPI_FLAVOR "MicrosoftHPC" CACHE STRING "Flavor of MPI." )

         # Use wmic to learn about the current machine
         execute_process(
            COMMAND wmic cpu get NumberOfCores
            OUTPUT_VARIABLE MPI_CORES_PER_CPU
            OUTPUT_STRIP_TRAILING_WHITESPACE )
        execute_process(
            COMMAND wmic computersystem get NumberOfLogicalProcessors
            OUTPUT_VARIABLE MPIEXEC_MAX_NUMPROCS
            OUTPUT_STRIP_TRAILING_WHITESPACE )
         execute_process(
            COMMAND wmic computersystem get NumberOfProcessors
            OUTPUT_VARIABLE MPI_CPUS_PER_NODE
            OUTPUT_STRIP_TRAILING_WHITESPACE )
         string( REGEX REPLACE ".*([0-9]+)" "\\1" MPI_CORES_PER_CPU ${MPI_CORES_PER_CPU})
         string( REGEX REPLACE ".*([0-9]+)" "\\1" MPIEXEC_MAX_NUMPROCS ${MPIEXEC_MAX_NUMPROCS})
         string( REGEX REPLACE ".*([0-9]+)" "\\1" MPI_CPUS_PER_NODE ${MPI_CPUS_PER_NODE})

         set( MPI_CPUS_PER_NODE ${MPI_CPUS_PER_NODE} CACHE STRING
            "Number of multi-core CPUs per node" FORCE )
         set( MPI_CORES_PER_CPU ${MPI_CORES_PER_CPU} CACHE STRING
            "Number of cores per cpu" FORCE )
         set( MPIEXEC_MAX_NUMPROCS ${MPIEXEC_MAX_NUMPROCS} CACHE STRING
            "Total number of available MPI ranks" FORCE )

         # Check for hyperthreading - This is important for reserving
         # threads for OpenMP tests...

         math( EXPR MPI_MAX_NUMPROCS_PHYSICAL "${MPI_CPUS_PER_NODE} * ${MPI_CORES_PER_CPU}" )
         if( "${MPI_MAX_NUMPROCS_PHYSICAL}" STREQUAL "${MPIEXEC_MAX_NUMPROCS}" )
            set( MPI_HYPERTHREADING "OFF" CACHE BOOL "Are we using hyperthreading?" FORCE )
         else()
            set( MPI_HYPERTHREADING "ON" CACHE BOOL "Are we using hyperthreading?" FORCE )
         endif()

         set( MPIEXEC_OMP_POSTFLAGS "-exitcodes"
            CACHE STRING "extra mpirun flags (list)." FORCE )
      endif()

    else()
      # Set DRACO_C4 and other variables
      setupDracoMPIVars()
   endif() # NOT "${DRACO_C4}" STREQUAL "SCALAR"

#   set( MPI_SETUP_DONE ON CACHE INTERNAL "Have we completed the MPI setup call?" )
   if( ${MPI_CXX_FOUND} )
      message(STATUS "Looking for MPI...${MPIEXEC}")
   else()
      message(STATUS "Looking for MPI...not found")
   endif()

   mark_as_advanced( MPI_FLAVOR MPIEXEC_OMP_POSTFLAGS MPI_LIBRARIES )

endmacro( setupMPILibrariesWindows )

#----------------------------------------------------------------------#
# End setupMPI.cmake
#----------------------------------------------------------------------#
