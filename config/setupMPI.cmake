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

#------------------------------------------------------------------------------#
# OpenMPI:
#------------------------------------------------------------------------------#
# Comments from Gabriel Rockefeller concerning OpenMPI 1.10.5 (email dated
# 2017 January 26):
#
# First, the multi-rank case:

# To run an MPI+threads calculation under OpenMPI 1.10.5, distributing multiple
# ranks across nodes (actually, sockets; more on this below) and using OpenMP
# threads to fill available cores, binding each thread to a specific core:

# % msub -I -lnodes=[number of nodes]
# % setenv OMP_PLACES threads
# % setenv NRANKS_PER_NODE [number of ranks per node]
# % setenv OMP_NUM_THREADS `echo ${SLURM_CPUS_ON_NODE} / ${NRANKS_PER_NODE} | bc`
# % mpirun -np `echo ${NRANKS_PER_NODE} \* ${SLURM_NNODES} | bc` --map-by socket:PE=${OMP_NUM_THREADS},span --bind-to core [application and other arguments ...]

# For example,

# % msub -I -lnodes=2
# % setenv OMP_PLACES threads
# % setenv NRANKS_PER_NODE 8
# % setenv OMP_NUM_THREADS `echo ${SLURM_CPUS_ON_NODE} / ${NRANKS_PER_NODE} | bc`
# % mpirun -np `echo ${NRANKS_PER_NODE} \* ${SLURM_NNODES} | bc` --map-by socket:PE=${OMP_NUM_THREADS},span --bind-to core ./xthi | sort -k 4,6
# Hello from rank 00, thread 00, on ml305.localdomain. (core affinity = 00)
# Hello from rank 00, thread 01, on ml305.localdomain. (core affinity = 01)
# Hello from rank 01, thread 00, on ml305.localdomain. (core affinity = 02)
# Hello from rank 01, thread 01, on ml305.localdomain. (core affinity = 03)
# Hello from rank 02, thread 00, on ml305.localdomain. (core affinity = 04)
# Hello from rank 02, thread 01, on ml305.localdomain. (core affinity = 05)
# Hello from rank 03, thread 00, on ml305.localdomain. (core affinity = 06)
# Hello from rank 03, thread 01, on ml305.localdomain. (core affinity = 07)
# Hello from rank 04, thread 00, on ml305.localdomain. (core affinity = 08)
# Hello from rank 04, thread 01, on ml305.localdomain. (core affinity = 09)
# Hello from rank 05, thread 00, on ml305.localdomain. (core affinity = 10)
# Hello from rank 05, thread 01, on ml305.localdomain. (core affinity = 11)
# Hello from rank 06, thread 00, on ml305.localdomain. (core affinity = 12)
# Hello from rank 06, thread 01, on ml305.localdomain. (core affinity = 13)
# Hello from rank 07, thread 00, on ml305.localdomain. (core affinity = 14)
# Hello from rank 07, thread 01, on ml305.localdomain. (core affinity = 15)
# Hello from rank 08, thread 00, on ml308.localdomain. (core affinity = 00)
# Hello from rank 08, thread 01, on ml308.localdomain. (core affinity = 01)
# Hello from rank 09, thread 00, on ml308.localdomain. (core affinity = 02)
# Hello from rank 09, thread 01, on ml308.localdomain. (core affinity = 03)
# Hello from rank 10, thread 00, on ml308.localdomain. (core affinity = 04)
# Hello from rank 10, thread 01, on ml308.localdomain. (core affinity = 05)
# Hello from rank 11, thread 00, on ml308.localdomain. (core affinity = 06)
# Hello from rank 11, thread 01, on ml308.localdomain. (core affinity = 07)
# Hello from rank 12, thread 00, on ml308.localdomain. (core affinity = 08)
# Hello from rank 12, thread 01, on ml308.localdomain. (core affinity = 09)
# Hello from rank 13, thread 00, on ml308.localdomain. (core affinity = 10)
# Hello from rank 13, thread 01, on ml308.localdomain. (core affinity = 11)
# Hello from rank 14, thread 00, on ml308.localdomain. (core affinity = 12)
# Hello from rank 14, thread 01, on ml308.localdomain. (core affinity = 13)
# Hello from rank 15, thread 00, on ml308.localdomain. (core affinity = 14)
# Hello from rank 15, thread 01, on ml308.localdomain. (core affinity = 15)

# For clarity, that mpirun command expands to

# % mpirun -np 16 --map-by socket:PE=2,span --bind-to core ./xthi | sort -k 4,6

# Note that you'll get an error if you try to run a combination of ranks and
# threads that don't fit evenly into sockets; probably the most common scenario
# that'll trigger this is running only one rank on a multi-socket node. There
# are other odd scenarios, like running an unusual number of ranks per node but
# a small number of threads per rank, where mapping by socket is probably the
# best thing to do, which is why I recommend it first.

# (As a digression, notice the affinity of threads on one node when I run two
# ranks, six threads per rank, and map by socket:

# % mpirun -np 2 --map-by socket:PE=6,span --bind-to core ./xthi | sort -k 4,6
# Hello from rank 00, thread 00, on ml305.localdomain. (core affinity = 00)
# Hello from rank 00, thread 01, on ml305.localdomain. (core affinity = 01)
# Hello from rank 00, thread 02, on ml305.localdomain. (core affinity = 02)
# Hello from rank 00, thread 03, on ml305.localdomain. (core affinity = 03)
# Hello from rank 00, thread 04, on ml305.localdomain. (core affinity = 04)
# Hello from rank 00, thread 05, on ml305.localdomain. (core affinity = 05)
# Hello from rank 01, thread 00, on ml305.localdomain. (core affinity = 08)
# Hello from rank 01, thread 01, on ml305.localdomain. (core affinity = 09)
# Hello from rank 01, thread 02, on ml305.localdomain. (core affinity = 10)
# Hello from rank 01, thread 03, on ml305.localdomain. (core affinity = 11)
# Hello from rank 01, thread 04, on ml305.localdomain. (core affinity = 12)
# Hello from rank 01, thread 05, on ml305.localdomain. (core affinity = 13)

# and compare the affinities when mapping by node instead:

# % mpirun -np 2 --map-by node:PE=6,span --bind-to core ./xthi | sort -k 4,6
# Hello from rank 00, thread 00, on ml305.localdomain. (core affinity = 00)
# Hello from rank 00, thread 01, on ml305.localdomain. (core affinity = 01)
# Hello from rank 00, thread 02, on ml305.localdomain. (core affinity = 02)
# Hello from rank 00, thread 03, on ml305.localdomain. (core affinity = 03)
# Hello from rank 00, thread 04, on ml305.localdomain. (core affinity = 04)
# Hello from rank 00, thread 05, on ml305.localdomain. (core affinity = 05)
# Hello from rank 01, thread 00, on ml305.localdomain. (core affinity = 06)
# Hello from rank 01, thread 01, on ml305.localdomain. (core affinity = 07)
# Hello from rank 01, thread 02, on ml305.localdomain. (core affinity = 08)
# Hello from rank 01, thread 03, on ml305.localdomain. (core affinity = 09)
# Hello from rank 01, thread 04, on ml305.localdomain. (core affinity = 10)
# Hello from rank 01, thread 05, on ml305.localdomain. (core affinity = 11)

# Cores 0-7 live on the first socket; cores 8-15 live on the second. Sockets
# represent "pools of memory bandwidth", so it's probably best to divide ranks
# (and their associated threads) evenly among sockets, using "--map-by socket",
# rather than filling the first socket and leaving the second underutilized, as
# happens when mapping by node.)

# To run a single-rank calculation and fill all available cores with OpenMP
# threads (as an example of a case where mapping by node is required):

# % msub -I -lnodes=1
# % setenv OMP_PLACES threads
# % setenv OMP_NUM_THREADS $SLURM_CPUS_ON_NODE
# % mpirun -np 1 --map-by node:PE=${OMP_NUM_THREADS},span --bind-to core ./xthi | sort -k 4,6
# Hello from rank 00, thread 00, on ml305.localdomain. (core affinity = 00)
# Hello from rank 00, thread 01, on ml305.localdomain. (core affinity = 01)
# Hello from rank 00, thread 02, on ml305.localdomain. (core affinity = 02)
# Hello from rank 00, thread 03, on ml305.localdomain. (core affinity = 03)
# Hello from rank 00, thread 04, on ml305.localdomain. (core affinity = 04)
# Hello from rank 00, thread 05, on ml305.localdomain. (core affinity = 05)
# Hello from rank 00, thread 06, on ml305.localdomain. (core affinity = 06)
# Hello from rank 00, thread 07, on ml305.localdomain. (core affinity = 07)
# Hello from rank 00, thread 08, on ml305.localdomain. (core affinity = 08)
# Hello from rank 00, thread 09, on ml305.localdomain. (core affinity = 09)
# Hello from rank 00, thread 10, on ml305.localdomain. (core affinity = 10)
# Hello from rank 00, thread 11, on ml305.localdomain. (core affinity = 11)
# Hello from rank 00, thread 12, on ml305.localdomain. (core affinity = 12)
# Hello from rank 00, thread 13, on ml305.localdomain. (core affinity = 13)
# Hello from rank 00, thread 14, on ml305.localdomain. (core affinity = 14)
# Hello from rank 00, thread 15, on ml305.localdomain. (core affinity = 15)

# Note the change from "--map-by socket" to "--map-by node".
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
#------------------------------------------------------------------------------#
# # Comments from Gabriel Rockefeller concerning OpenMPI 1.10.5 (email dated
# 2017 January 26):
#
# For a MPI+OpenMP job that has exclusive access to its node:
#
# % setenv OMP_NUM_THREADS [number of threads per rank]
# % setenv OMP_PLACES threads
# % aprun -n [num_ranks] -N [ranks_per_node] --cc depth -d $OMP_NUM_THREADS -j 1 [your binary plus options]
#
# where ranks_per_node * OMP_NUM_THREADS should probably come out to the number
# of PEs on each node (unless you're intentionally undersubscribing each node).
#
# (Optionally, if you are undersubscribing the node, add -S
# [ranks_per_numa_node] to spread ranks out across both NUMA nodes; whether this
# helps or hurts will depend on details of the particular calculation you're
# running, i.e., whether communication between ranks or access to local memory
# is more important. if you add it, ranks_per_numa_node should probably be
# ranks_per_node / 2, on Trinitite/Trinity Haswell.)
#
# For example, on Trinitite:
#
# % msub -I -lnodes=1,walltime=1:00:00
# % setenv OMP_NUM_THREADS 4
# % setenv OMP_PLACES threads
# % module load xthi
# % aprun -n 4 --cc depth -d $OMP_NUM_THREADS -j 1 xthi | sort -k 4,6
# Hello from rank 00, thread 00, on nid00030. (core affinity = 00)
# Hello from rank 00, thread 01, on nid00030. (core affinity = 01)
# Hello from rank 00, thread 02, on nid00030. (core affinity = 02)
# Hello from rank 00, thread 03, on nid00030. (core affinity = 03)
# Hello from rank 01, thread 00, on nid00030. (core affinity = 04)
# Hello from rank 01, thread 01, on nid00030. (core affinity = 05)
# Hello from rank 01, thread 02, on nid00030. (core affinity = 06)
# Hello from rank 01, thread 03, on nid00030. (core affinity = 07)
# Hello from rank 02, thread 00, on nid00030. (core affinity = 08)
# Hello from rank 02, thread 01, on nid00030. (core affinity = 09)
# Hello from rank 02, thread 02, on nid00030. (core affinity = 10)
# Hello from rank 02, thread 03, on nid00030. (core affinity = 11)
# Hello from rank 03, thread 00, on nid00030. (core affinity = 12)
# Hello from rank 03, thread 01, on nid00030. (core affinity = 13)
# Hello from rank 03, thread 02, on nid00030. (core affinity = 14)
# Hello from rank 03, thread 03, on nid00030. (core affinity = 15)
#
# The same, with -S 2, puts two ranks on each NUMA node, instead of all four
# ranks on the first NUMA node:
#
# % aprun -n 4 -S 2 --cc depth -d $OMP_NUM_THREADS -j 1 xthi | sort -k 4,6
# Hello from rank 00, thread 00, on nid00030. (core affinity = 00)
# Hello from rank 00, thread 01, on nid00030. (core affinity = 01)
# Hello from rank 00, thread 02, on nid00030. (core affinity = 02)
# Hello from rank 00, thread 03, on nid00030. (core affinity = 03)
# Hello from rank 01, thread 00, on nid00030. (core affinity = 04)
# Hello from rank 01, thread 01, on nid00030. (core affinity = 05)
# Hello from rank 01, thread 02, on nid00030. (core affinity = 06)
# Hello from rank 01, thread 03, on nid00030. (core affinity = 07)
# Hello from rank 02, thread 00, on nid00030. (core affinity = 16)
# Hello from rank 02, thread 01, on nid00030. (core affinity = 17)
# Hello from rank 02, thread 02, on nid00030. (core affinity = 18)
# Hello from rank 02, thread 03, on nid00030. (core affinity = 19)
# Hello from rank 03, thread 00, on nid00030. (core affinity = 20)
# Hello from rank 03, thread 01, on nid00030. (core affinity = 21)
# Hello from rank 03, thread 02, on nid00030. (core affinity = 22)
# Hello from rank 03, thread 03, on nid00030. (core affinity = 23)
#
# In contrast, for multiple jobs on the same node, unset OMP_PLACES and use --cc
# none to avoid binding threads to particular PEs or ranges:
#
# % unsetenv OMP_PLACES
# % aprun -n 4 --cc none -d $OMP_NUM_THREADS -j 1 xthi | sort -k 4,6
# Hello from rank 00, thread 00, on nid00030. (core affinity = 00-63)
# Hello from rank 00, thread 01, on nid00030. (core affinity = 00-63)
# Hello from rank 00, thread 02, on nid00030. (core affinity = 00-63)
# Hello from rank 00, thread 03, on nid00030. (core affinity = 00-63)
# Hello from rank 01, thread 00, on nid00030. (core affinity = 00-63)
# Hello from rank 01, thread 01, on nid00030. (core affinity = 00-63)
# Hello from rank 01, thread 02, on nid00030. (core affinity = 00-63)
# Hello from rank 01, thread 03, on nid00030. (core affinity = 00-63)
# Hello from rank 02, thread 00, on nid00030. (core affinity = 00-63)
# Hello from rank 02, thread 01, on nid00030. (core affinity = 00-63)
# Hello from rank 02, thread 02, on nid00030. (core affinity = 00-63)
# Hello from rank 02, thread 03, on nid00030. (core affinity = 00-63)
# Hello from rank 03, thread 00, on nid00030. (core affinity = 00-63)
# Hello from rank 03, thread 01, on nid00030. (core affinity = 00-63)
# Hello from rank 03, thread 02, on nid00030. (core affinity = 00-63)
# Hello from rank 03, thread 03, on nid00030. (core affinity = 00-63)
#
# I'll add this to Confluence once you've had a chance to try it. Running
# multiple jobs on a Trinity/Trinitite node involves Adaptive's MAPN
# functionality, which might not quite be working as intended yet.
#
#------------------------------------------------------------------------------#
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

  # -b        Bypass transfer of application executable to the compute node.
  # -cc none  Do not bind threads to a CPU within the assigned NUMA node.
  # -q        Quiet
  # -m 1400m  Reserve 1.4 GB of RAM per PE. Trinitite/Trinity has 4GB/core for
  #           haswells, 1.4GB/core for KNL
  # -F shared enabled shared mode to allow multiple applications to run on a
  #           single node.
  # set( MPIEXEC_POSTFLAGS "-q -F shared -b -m 1400m" CACHE STRING
  #   "extra mpirun flags (list)." FORCE)
   set( MPIEXEC_POSTFLAGS "-N 1 --cpu_bind=cores --gres=craynetwork:0 --exclusive"
     CACHE STRING
     "extra mpirun flags (list)." FORCE)
    # Consider using 'aprun -n # -N # -S # -d # -T -cc depth ...'
    # -n #  number of processes
    # -N #  number of processes per node
    # -S #  number of processes per numa node
    # -d #  cpus-per-pe
    # -T    sync-output
    # -cc depth PEs are constrained to CPUs with a distance of depth between
    #       them so each PE's threads can be constrained to the CPUs closest to
    #       the PE's CPU.
    #set( MPIEXEC_OMP_POSTFLAGS "-q -b -d $ENV{OMP_NUM_THREADS}"
    #  CACHE STRING "extra mpirun flags (list)." FORCE)

    set( MPIEXEC_OMP_POSTFLAGS "${MPIEXEC_POSTFLAGS} -c ${MPI_CORES_PER_CPU}"
      CACHE STRING "extra mpirun flags (list)." FORCE)
  # Extra flags for OpenMP + MPI
#   if( DEFINED ENV{OMP_NUM_THREADS} )

#   else()
#     message( STATUS "
# WARNING: ENV{OMP_NUM_THREADS} is not set in your environment,
#          all OMP tests will be disabled." )
#   endif()

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

      # elseif( "${MPIEXEC}" MATCHES aprun)
      elseif( CRAY_PE )
        setupCrayMPI()

      elseif( "${MPIEXEC}" MATCHES srun)
        setupSequoiaMPI()

      else()
         message( FATAL_ERROR "
The Draco build system doesn't know how to configure the build for
  MPIEXEC     = ${MPIEXEC}
  DBS_MPI_VER = ${DBS_MPI_VER}")
      endif()

      # Mark some of the variables created by the above logic as
      # 'advanced' so that they do not show up in the 'simple' ccmake
      # view.
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
