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
         if( NOT MPIEXEC )
            if( "${SITE}" MATCHES "c[it]" ) 
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
      get_filename_component( compiler_wo_path "${CMAKE_CXX_COMPILER}" NAME )
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
               
         # Ref: http://www.open-mpi.org/faq/?category=tuning#using-paffinity-v1.2
         # The --bind-to-none option available (and the default) in
         # OpenMPI 1.4+ is apparently not a synonym for this;
         #
         # % ompi_info -param mpi all
         #
         # with OpenMPI 1.6.3 on Moonlight reports that
         # mpi_paffinity_alone is 1, which means, "[A]ssume that this
         # job is the only (set of) process(es) running on each node
         # and bind processes to processors, starting with processor
         # ID 0."  Setting mpi_paffinity_alone to 0 allows parallel
         # ctest to work correctly.  MPIEXEC_POSTFLAGS only affects
         # MPI-only tests (and not MPI+OpenMP tests).

         # This flag also shows up in
         # jayenne/pkg_tools/run_milagro_test.py and regress_funcs.py.
         if( ${DBS_MPI_VER_MAJOR}.${DBS_MPI_VER_MINOR} VERSION_LESS 1.7 )
           set( MPIEXEC_POSTFLAGS --mca mpi_paffinity_alone 0 )
           set( MPIEXEC_POSTFLAGS_STRING "--mca mpi_paffinity_alone 0" )
         else()
           # Flag provided by Sam Gutierrez (2014-04-08),
           set( MPIEXEC_POSTFLAGS -mca hwloc_base_binding_policy none )
           set( MPIEXEC_POSTFLAGS_STRING "-mca hwloc_base_binding_policy none" )
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
         math( EXPR MPI_CPUS_PER_NODE 
           "${MPIEXEC_MAX_NUMPROCS} / ${MPI_CORES_PER_CPU}" )
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
            set( MPI_HYPERTHREADING "OFF" CACHE BOOL 
              "Are we using hyperthreading?" FORCE )
         else()
            set( MPI_HYPERTHREADING "ON" CACHE BOOL 
              "Are we using hyperthreading?" FORCE )
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
           #message( "MPI Version < 1.7" )
           set( MPIEXEC_OMP_POSTFLAGS 
             -bind-to-socket -cpus-per-proc ${MPI_CORES_PER_CPU} --report-bindings
             CACHE STRING "extra mpirun flags (list)." FORCE )
           set( MPIEXEC_OMP_POSTFLAGS_STRING 
             "-bind-to-socket -cpus-per-proc ${MPI_CORES_PER_CPU} --report-bindings"
             CACHE STRING "extra mpirun flags (list)." FORCE)
         elseif( ${DBS_MPI_VER_MAJOR}.${DBS_MPI_VER_MINOR} VERSION_GREATER 1.7 )
           #message( "MPI Version > 1.7" )
           set( MPIEXEC_OMP_POSTFLAGS 
             -bind-to socket --map-by ppr:${MPI_CORES_PER_CPU}:socket --report-bindings
             CACHE STRING "extra mpirun flags (list)." FORCE )
           set( MPIEXEC_OMP_POSTFLAGS_STRING 
             "-bind-to socket --map-by ppr:${MPI_CORES_PER_CPU}:socket --report-bindings"
             CACHE STRING "extra mpirun flags (list)." FORCE)
        else() 
          #message( "MPI Version = 1.7" )
          # Version 1.7.4
           set( MPIEXEC_OMP_POSTFLAGS 
             -bind-to socket --map-by socket:PPR=${MPI_CORES_PER_CPU}
             CACHE STRING "extra mpirun flags (list)." FORCE )
           set( MPIEXEC_OMP_POSTFLAGS_STRING 
             "-bind-to socket --map-by socket:PPR=${MPI_CORES_PER_CPU}"
             CACHE STRING "extra mpirun flags (list)." FORCE)
           
        endif()
         mark_as_advanced( MPI_FLAVOR MPIEXEC_OMP_POSTFLAGS_STRING MPIEXEC_OMP_POSTFLAGS
            MPI_LIBRARIES )         

      ### Cray wrappers for mpirun
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

      ### Sequoia
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

      ### Intel MPI
      elseif( "${MPIEXEC}" MATCHES intel-mpi OR
          "${DBS_MPI_VER}" MATCHES "Intel[(]R[)] MPI Library" )

         set( MPI_FLAVOR "impi" CACHE STRING "Flavor of MPI." )

         # Find the version of Intel MPI

         if( "${DBS_MPI_VER}" MATCHES "[0-9].[0-9].[0-9]" )
            string( REGEX REPLACE ".*([0-9])[.]([0-9])[.]([0-9]).*" "\\1"
               DBS_MPI_VER_MAJOR ${DBS_MPI_VER} )
            string( REGEX REPLACE ".*([0-9])[.]([0-9])[.]([0-9]).*" "\\2"
               DBS_MPI_VER_MINOR ${DBS_MPI_VER} )
         elseif( "${DBS_MPI_VER}" MATCHES "[0-9].[0-9]" )
            string( REGEX REPLACE ".*([0-9]).([0-9]).*" "\\1"
               DBS_MPI_VER_MAJOR ${DBS_MPI_VER} )
            string( REGEX REPLACE ".*([0-9]).([0-9]).*" "\\2"
               DBS_MPI_VER_MINOR ${DBS_MPI_VER} )
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
         math( EXPR MPI_CPUS_PER_NODE 
           "${MPIEXEC_MAX_NUMPROCS} / ${MPI_CORES_PER_CPU}" )
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
            set( MPI_HYPERTHREADING "OFF" CACHE BOOL 
              "Are we using hyperthreading?" FORCE )
         else()
            set( MPI_HYPERTHREADING "ON" CACHE BOOL 
              "Are we using hyperthreading?" FORCE )
         endif()

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
      # (version ${DBS_MPI_VER})")

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

##---------------------------------------------------------------------------##
## 
##---------------------------------------------------------------------------##

macro( setupMPILibrariesWindows )

   # MPI ---------------------------------------------------------------------
   if( NOT "${DRACO_C4}" STREQUAL "SCALAR" ) # AND "${MPIEXEC}x" STREQUAL "x" )

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
             -bind-to-socket -cpus-per-proc ${MPI_CORES_PER_CPU} --report-bindings
             CACHE STRING "extra mpirun flags (list)." FORCE )
           set( MPIEXEC_OMP_POSTFLAGS_STRING 
             "-bind-to-socket -cpus-per-proc ${MPI_CORES_PER_CPU} --report-bindings"
             CACHE STRING "extra mpirun flags (list)." FORCE)
         else()
           set( MPIEXEC_OMP_POSTFLAGS 
             -bind-to socket --map-by ppr:${MPI_CORES_PER_CPU}:socket --report-bindings
             CACHE STRING "extra mpirun flags (list)." FORCE )
           set( MPIEXEC_OMP_POSTFLAGS_STRING 
             "-bind-to socket --map-by ppr:${MPI_CORES_PER_CPU}:socket --report-bindings"
             CACHE STRING "extra mpirun flags (list)." FORCE)
         endif()
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
         set( MPIEXEC_OMP_POSTFLAGS -exitcodes            
            CACHE STRING "extra mpirun flags (list)." FORCE )
         set( MPIEXEC_OMP_POSTFLAGS_STRING "-exitcodes"
            CACHE STRING "extra mpirun flags (list)." FORCE)
         mark_as_advanced( MPI_FLAVOR MPIEXEC_OMP_POSTFLAGS_STRING MPIEXEC_OMP_POSTFLAGS
            MPI_LIBRARIES )    
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
         set( MPIEXEC_OMP_POSTFLAGS -exitcodes            
            CACHE STRING "extra mpirun flags (list)." FORCE )
         set( MPIEXEC_OMP_POSTFLAGS_STRING "-exitcodes"
            CACHE STRING "extra mpirun flags (list)." FORCE)
         mark_as_advanced( MPI_FLAVOR MPIEXEC_OMP_POSTFLAGS_STRING MPIEXEC_OMP_POSTFLAGS
            MPI_LIBRARIES )                
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
