#-----------------------------*-cmake-*----------------------------------------#
# file   config/setupMPI.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   2016 Sep 22
# brief  Setup MPI Vendors
# note   Copyright (C) 2016-2019 Triad National Security, LLC.
#        All rights reserved.
#
# Try to find MPI in the default locations (look for mpic++ in PATH)
#
# See cmake --help-module FindMPI for details on variables set and published
# targets. Additionally, this module will set the following variables:
#
# DRACO_C4   MPI|SCALAR
# C4_SCALAR  BOOL
# C4_MPI     BOOL
#
#------------------------------------------------------------------------------#

include_guard(GLOBAL)
include( FeatureSummary )

##---------------------------------------------------------------------------##
## Set MPI flavor and vendor version
##
## Returns (as cache variables)
## - MPI_VERSION
## - MPI_FLAVOR = {openmpi, mpich, cray, spectrum, mvapich2, intel}
##---------------------------------------------------------------------------##
function( setMPIflavorVer )

  # First attempt to determine MPI flavor -- scape flavor from full path
  # (this ususally works for HPC or systems with modules)
  if( CRAY_PE )
    set( MPI_FLAVOR "cray" )
  elseif( "${MPIEXEC_EXECUTABLE}" MATCHES "openmpi")
    set( MPI_FLAVOR "openmpi" )
  elseif( "${MPIEXEC_EXECUTABLE}" MATCHES "mpich")
    set( MPI_FLAVOR "mpich" )
  elseif( "${MPIEXEC_EXECUTABLE}" MATCHES "impi" OR
      "${MPIEXEC_EXECUTABLE}" MATCHES "clusterstudio" )
    set( MPI_FLAVOR "intel" )
  elseif( "${MPIEXEC_EXECUTABLE}" MATCHES "mvapich2")
    set( MPI_FLAVOR "mvapich2" )
  elseif( "${MPIEXEC_EXECUTABLE}" MATCHES "spectrum-mpi" OR
      "${MPIEXEC_EXECUTABLE}" MATCHES "jsrun" )
    set( MPI_FLAVOR "spectrum")
  endif()

  if( CRAY_PE )
    if( DEFINED ENV{CRAY_MPICH2_VER} )
      set( MPI_VERSION $ENV{CRAY_MPICH2_VER} )
    endif()
  else()
    execute_process( COMMAND ${MPIEXEC_EXECUTABLE} --version
      OUTPUT_VARIABLE DBS_MPI_VER_OUT
      ERROR_VARIABLE DBS_MPI_VER_ERR)
    set( DBS_MPI_VER "${DBS_MPI_VER_OUT}${DBS_MPI_VER_ERR}")
    string(REPLACE "\n" ";" TEMP ${DBS_MPI_VER})
    foreach( line ${TEMP} )

      # extract the version...
      if( ${line} MATCHES "Version" OR ${line} MATCHES "OpenRTE" OR
          ${line} MATCHES "Open MPI" OR ${line} MATCHES "Spectrum MPI")
        set(DBS_MPI_VER "${line}")
        if( "${DBS_MPI_VER}" MATCHES "[0-9]+[.][0-9]+[.][0-9]+" )
          string( REGEX REPLACE ".*([0-9]+)[.]([0-9]+)[.]([0-9]+).*" "\\1"
            DBS_MPI_VER_MAJOR ${DBS_MPI_VER} )
          string( REGEX REPLACE ".*([0-9]+)[.]([0-9]+)[.]([0-9]+).*" "\\2"
            DBS_MPI_VER_MINOR ${DBS_MPI_VER} )
          string( REGEX REPLACE ".*([0-9]+)[.]([0-9]+)[.]([0-9]+).*" "\\3"
            DBS_MPI_VER_PATCH ${DBS_MPI_VER} )
          set( MPI_VERSION
            "${DBS_MPI_VER_MAJOR}.${DBS_MPI_VER_MINOR}.${DBS_MPI_VER_PATCH}" )
        elseif( "${DBS_MPI_VER}" MATCHES "[0-9]+[.][0-9]+" )
          string( REGEX REPLACE ".*([0-9]+)[.]([0-9]+).*" "\\1"
            DBS_MPI_VER_MAJOR ${DBS_MPI_VER} )
          string( REGEX REPLACE ".*([0-9]+)[.]([0-9]+).*" "\\2"
            DBS_MPI_VER_MINOR ${DBS_MPI_VER} )
          set( MPI_VERSION "${DBS_MPI_VER_MAJOR}.${DBS_MPI_VER_MINOR}" )
        endif()
      else()
        message(FATAL_ERROR "DBS did not find the MPI version string (is this "
          "an older openmpi?)")
      endif()

      # if needed, make a 2nd pass at identifying the MPI flavor
      if( NOT DEFINED MPI_FLAVOR )
        if ("${line}" MATCHES "HYDRA")
          set(MPI_FLAVOR "mpich")
        elseif("${line}" MATCHES "OpenRTE")
          set(MPI_FLAVOR "openmpi")
        elseif("${line}" MATCHES "intel-mpi" OR
            "${line}" MATCHES "Intel[(]R[)] MPI Library" )
          set(MPI_FLAVOR "intel")
        endif()
      endif()

      # Once we have the needed information stop parsing...
      if( DEFINED MPI_FLAVOR AND DEFINED MPI_VERSION )
        break()
      endif()
    endforeach()

  endif()

  # if the FindMPI.cmake module didn't set the version, then try to do so here.
  if( NOT DEFINED MPI_VERSION AND DEFINED MPI_C_VERSION)
    set( MPI_VERSION ${MPI_C_VERSION} )
  endif()

  # Return the discovered values to the calling scope
  if( DEFINED MPI_FLAVOR )
    set( MPI_FLAVOR "${MPI_FLAVOR}" CACHE STRING "Vendor brand of MPI" FORCE )
  endif()
  if( DEFINED MPI_VERSION )
    set( MPI_VERSION "${MPI_VERSION}" CACHE STRING "Vendor version of MPI"
      FORCE )
  endif()

endfunction()

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

  # These cmake commands, while useful, don't provide the topology detail that we
  # are interested in (i.e. number of sockets per node). We could use the results
  # of these queries to know if hyperthreading is enabled (if logical != physical
  # cores)
  # - cmake_host_system_information(RESULT MPI_PHYSICAL_CORES
  #   QUERY NUMBER_OF_PHYSICAL_CORES)
  # - cmake_host_system_information(RESULT MPI_LOGICAL_CORES
  #   QUERY NUMBER_OF_LOGICAL_CORES)

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
  math( EXPR MPI_MAX_NUMPROCS_PHYSICAL
    "${MPI_PHYSICAL_CORES} * ${MPI_CORES_PER_CPU}" )
  if( "${MPI_MAX_NUMPROCS_PHYSICAL}" STREQUAL "${MPIEXEC_MAX_NUMPROCS}" )
    set( MPI_HYPERTHREADING "OFF" CACHE BOOL "Are we using hyperthreading?"
      FORCE )
  else()
    set( MPI_HYPERTHREADING "ON" CACHE BOOL "Are we using hyperthreading?"
      FORCE )
  endif()
endmacro()

##---------------------------------------------------------------------------##
## Setup OpenMPI
##---------------------------------------------------------------------------##
macro( setupOpenMPI )

  # sanity check, these OpenMPI flags (below) require version >= 1.8
  if( ${MPI_VERSION} VERSION_LESS 1.8 )
    message( FATAL_ERROR "OpenMPI version < 1.8 found." )
  endif()

  # Setting mpi_paffinity_alone to 0 allows parallel ctest to work correctly.
  # MPIEXEC_POSTFLAGS only affects MPI-only tests (and not MPI+OpenMP tests).
  # . --oversubscribe is only available for openmpi version >= 3.0
  # . -H localhost,localhost,localhost,localhost might work for older versions.
  if( "$ENV{GITLAB_CI}" STREQUAL "true" OR "$ENV{TRAVIS}" STREQUAL "true")
    set(runasroot "--allow-run-as-root --oversubscribe")
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
    # "-bind-to socket --map-by ppr:${MPI_CORES_PER_CPU}:socket
    # --report-bindings ${runasroot}"
  endif()

  set( MPIEXEC_OMP_POSTFLAGS ${MPIEXEC_OMP_POSTFLAGS}
    CACHE STRING "extra mpirun flags (list)." FORCE )

  mark_as_advanced( MPI_CPUS_PER_NODE MPI_CORES_PER_CPU
    MPI_PHYSICAL_CORES MPI_MAX_NUMPROCS_PHYSICAL MPI_HYPERTHREADING )

endmacro()

##---------------------------------------------------------------------------##
## Setup Mpich
##---------------------------------------------------------------------------##
macro( setupMpichMPI )

  # Find cores/cpu, cpu/node, hyperthreading
  query_topology()

endmacro()

##---------------------------------------------------------------------------##
## Setup Intel MPI
##---------------------------------------------------------------------------##
macro( setupIntelMPI )

  # Find cores/cpu, cpu/node, hyperthreading
  query_topology()

endmacro()

##---------------------------------------------------------------------------##
## Setup Cray MPI wrappers (APRUN)
##---------------------------------------------------------------------------##
macro( setupCrayMPI )

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
  string(APPEND postflags " --exclusive --gres=craynetwork:0")
  set( MPIEXEC_POSTFLAGS ${postflags} CACHE STRING
    "extra mpirun flags (list)." FORCE)

  set( MPIEXEC_OMP_POSTFLAGS "${MPIEXEC_POSTFLAGS} -c ${MPI_CORES_PER_CPU}"
    CACHE STRING "extra mpirun flags (list)." FORCE)

endmacro()

##---------------------------------------------------------------------------##
## Setup Spectrum MPI wrappers (Sierra, Rzansel, Rzmanta, Ray)
##---------------------------------------------------------------------------##
macro( setupSpectrumMPI )

  # Find cores/cpu, cpu/node, hyperthreading
  query_topology()

endmacro()

#------------------------------------------------------------------------------#
# Setup MPI when on Linux
#------------------------------------------------------------------------------#
macro( setupMPILibrariesUnix )

  # MPI ---------------------------------------------------------------------
  if( NOT "${DRACO_C4}" STREQUAL "SCALAR" )

    message(STATUS "Looking for MPI...")

    # Preserve data that may already be set.
    if( DEFINED ENV{MPIRUN} )
      set( MPIEXEC_EXECUTABLE $ENV{MPIRUN} CACHE STRING
        "Program to execute MPI parallel programs." )
    elseif( DEFINED ENV{MPIEXEC_EXECUTABLE} )
      set( MPIEXEC_EXECUTABLE $ENV{MPIEXEC_EXECUTABLE} CACHE STRING
        "Program to execute MPI parallel programs." )
    elseif( DEFINED ENV{MPIEXEC} )
      set( MPIEXEC_EXECUTABLE $ENV{MPIEXEC} CACHE STRING
        "Program to execute MPI parallel programs." )
    endif()

    # If this is a Cray system and the Cray MPI compile wrappers are used,
    # then do some special setup:

    if( CRAY_PE )
      if( NOT EXISTS ${MPIEXEC_EXECUTABLE} )
        find_program( MPIEXEC_EXECUTABLE srun )
      endif()
      set( MPIEXEC_EXECUTABLE ${MPIEXEC_EXECUTABLE} CACHE STRING
        "Program to execute MPI parallel programs." FORCE )
      set( MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING
        "mpirun flag used to specify the number of processors to use")
    endif()

    # Call the standard CMake FindMPI macro.
    find_package( MPI QUIET )

    # Try to discover the MPI flavor and the vendor version. Returns
    # MPI_VERSION, MPI_FLAVOR as cache variables
    setMPIflavorVer()

    # Set additional flags, environments that are MPI vendor specific.
    if( "${MPI_FLAVOR}" MATCHES "openmpi" )
      setupOpenMPI()
    elseif( "${MPI_FLAVOR}" MATCHES "mpich" )
      setupMpichMPI()
    elseif( "${MPI_FLAVOR}" MATCHES "intel" )
      setupIntelMPI()
    elseif( "${MPI_FLAVOR}" MATCHES "spectrum" )
      setupSpectrumMPI()
    elseif( "${MPI_FLAVOR}" MATCHES "cray" )
      setupCrayMPI()
    else()
      message( FATAL_ERROR "
The Draco build system doesn't know how to configure the build for
  MPIEXEC_EXECUTABLE     = ${MPIEXEC_EXECUTABLE}
  DBS_MPI_VER = ${DBS_MPI_VER}
  CRAY_PE     = ${CRAY_PE}")
    endif()

    # Mark some of the variables created by the above logic as 'advanced' so
    # that they do not show up in the 'simple' ccmake view.
    mark_as_advanced( MPI_EXTRA_LIBRARY MPI_LIBRARY )

    message(STATUS "Looking for MPI.......found ${MPIEXEC_EXECUTABLE}")

    # Sanity Checks for DRACO_C4==MPI
    if( "${MPI_CORES_PER_CPU}x" STREQUAL "x" )
      message( FATAL_ERROR "setupMPILibrariesUnix:: MPI_CORES_PER_CPU "
        "is not set!")
    endif()

  endif()

  # Set DRACO_C4 and other variables. Returns DRACO_C4, C4_SCALAR, C4_MPI
  setupDracoMPIVars()

  set_package_properties( MPI PROPERTIES
    URL "http://www.open-mpi.org/"
    DESCRIPTION "A High Performance Message Passing Library"
    TYPE RECOMMENDED
    PURPOSE "If not available, all Draco components will be built as scalar applications."  )

  mark_as_advanced( MPIEXEC_OMP_POSTFLAGS MPI_LIBRARIES )

endmacro()

##---------------------------------------------------------------------------##
## setupMPILibrariesWindows
##---------------------------------------------------------------------------##

macro( setupMPILibrariesWindows )

  set( verbose FALSE )

  # MPI ---------------------------------------------------------------------
  if( NOT "${DRACO_C4}" STREQUAL "SCALAR" )

    message(STATUS "Looking for MPI...")
    find_package( MPI QUIET )

    # If this macro is called from a MinGW builds system (for a CAFS
    # subdirectory) and is trying to MS-MPI, the above check will fail (when
    # cmake > 3.12). However, MS-MPI is known to be good when linking with
    # Visual Studio so override the 'failed' report.
    if( "${MPI_C_LIBRARIES}" MATCHES "msmpi" AND
        "${CMAKE_GENERATOR}" STREQUAL "MinGW Makefiles")
      if( EXISTS "${MPI_C_LIBRARIES}" AND EXISTS "${MPI_C_INCLUDE_DIRS}" )
        set( MPI_C_FOUND TRUE )
        set( MPI_Fortran_FOUND TRUE )
      endif()
    endif()

    # For MS-MPI, mpifptr.h is architecture dependent. Figure out what arch this
    # is and save this path to MPI_Fortran_INCLUDE_PATH
    list( GET MPI_C_LIBRARIES 0 first_c_mpi_library )
    if( first_c_mpi_library AND NOT MPI_Fortran_INCLUDE_PATH )
      get_filename_component( MPI_Fortran_INCLUDE_PATH
        "${first_c_mpi_library}" DIRECTORY )
      string( REGEX REPLACE "[Ll]ib" "Include" MPI_Fortran_INCLUDE_PATH
        ${MPI_Fortran_INCLUDE_PATH} )
      set( MPI_Fortran_INCLUDE_PATH
        "${MPI_CXX_INCLUDE_PATH};${MPI_Fortran_INCLUDE_PATH}"
        CACHE STRING "Location for MPI include files for Fortran.")
      if( verbose )
        message("MPI_Fortran_INCLUDE_PATH = ${MPI_Fortran_INCLUDE_PATH}")
      endif()
    endif()

    setupDracoMPIVars()

    # Find the version
    # This is not working (hardwire it for now)
    execute_process( COMMAND "${MPIEXEC_EXECUTABLE}" -help
      OUTPUT_VARIABLE DBS_MPI_VER_OUT
      ERROR_VARIABLE DBS_MPI_VER_ERR
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_STRIP_TRAILING_WHITESPACE
      )
    if( "${DBS_MPI_VER_OUT}" MATCHES "Microsoft MPI Startup Program" )
      string( REGEX REPLACE ".*Version ([0-9.]+).*" "\\1" DBS_MPI_VER
        "${DBS_MPI_VER_OUT}${DBS_MPI_VER_ERR}")
      string( REGEX REPLACE ".*([0-9])[.]([0-9])[.]([0-9]+).*" "\\1"
        DBS_MPI_VER_MAJOR ${DBS_MPI_VER} )
      string( REGEX REPLACE ".*([0-9])[.]([0-9])[.]([0-9]+).*" "\\2"
        DBS_MPI_VER_MINOR ${DBS_MPI_VER} )
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
    if("${MPIEXEC_EXECUTABLE}" MATCHES "Microsoft MPI")
      set( MPI_FLAVOR "MicrosoftMPI" CACHE STRING "Flavor of MPI." )

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
      string( REGEX REPLACE ".*[\n]([0-9]+$)" "\\1" MPI_CORES_PER_CPU
        ${MPI_CORES_PER_CPU})
      string( REGEX REPLACE ".*[\n]([0-9]+$)" "\\1" MPIEXEC_MAX_NUMPROCS
        ${MPIEXEC_MAX_NUMPROCS})
      string( REGEX REPLACE ".*[\n]([0-9]+$)" "\\1" MPI_CPUS_PER_NODE
        ${MPI_CPUS_PER_NODE})

      set( MPI_CPUS_PER_NODE ${MPI_CPUS_PER_NODE} CACHE STRING
        "Number of multi-core CPUs per node" FORCE )
      set( MPI_CORES_PER_CPU ${MPI_CORES_PER_CPU} CACHE STRING
        "Number of cores per cpu" FORCE )
      set( MPIEXEC_MAX_NUMPROCS ${MPIEXEC_MAX_NUMPROCS} CACHE STRING
        "Total number of available MPI ranks" FORCE )

      # Check for hyperthreading - This is important for reserving
      # threads for OpenMP tests...

      math( EXPR MPI_MAX_NUMPROCS_PHYSICAL
        "${MPI_CPUS_PER_NODE} * ${MPI_CORES_PER_CPU}" )
      if( "${MPI_MAX_NUMPROCS_PHYSICAL}" STREQUAL "${MPIEXEC_MAX_NUMPROCS}" )
        set( MPI_HYPERTHREADING "OFF" CACHE BOOL
          "Are we using hyperthreading?" FORCE )
      else()
        set( MPI_HYPERTHREADING "ON" CACHE BOOL
          "Are we using hyperthreading?" FORCE )
      endif()

      set( MPIEXEC_OMP_POSTFLAGS "-exitcodes"
        CACHE STRING "extra mpirun flags (list)." FORCE )
    endif()

  else()
    # Set DRACO_C4 and other variables
    setupDracoMPIVars()
  endif() # NOT "${DRACO_C4}" STREQUAL "SCALAR"

  # Found MPI_C, but not MPI_CXX -- create a duplicate to satisfy link targets.
  if( TARGET MPI::MPI_C AND NOT TARGET MPI::MPI_CXX)

    # Windows systems with dll libraries.
    add_library( MPI::MPI_CXX SHARED IMPORTED )

    # Windows with dlls, but only Release libraries.
    set_target_properties(MPI::MPI_CXX PROPERTIES
      IMPORTED_LOCATION_RELEASE         "${MPI_C_LIBRARIES}"
      IMPORTED_IMPLIB                   "${MPI_C_LIBRARIES}"
      INTERFACE_INCLUDE_DIRECTORIES     "${MPI_C_INCLUDE_DIRS}"
      IMPORTED_CONFIGURATIONS           Release
      IMPORTED_LINK_INTERFACE_LANGUAGES "CXX" )

  endif()

  # Don't link to the C++ MS-MPI library when compiling with MinGW gfortran.
  # Instead, link to libmsmpi.a that was created via gendef.exe and dlltool.exe
  # from msmpi.dll.  Ref: http://www.geuz.org/pipermail/getdp/2012/001520.html,
  # or
  # https://github.com/KineticTheory/Linux-HPC-Env/wiki/Setup-Win32-development-environment

  # Preparing Microsoft's MPI to work with x86_64-w64-mingw32-gfortran by
  # creating libmsmpi.a. (Last tested: 2017-08-31)
  #
  # 1.) You need MinGW/MSYS. Please make sure that the Devel tools are
  #     installed.
  # 2.) Download and install Microsoft's MPI. You need the main program and
  #     the SDK.
  # 3.) In the file %MSMPI_INC%\mpif.h, replace INT_PTR_KIND() by 8
  # 4.) Create a MSYS version of the MPI library:
  #     cd %TEMP%
  #     copy c:\Windows\System32\msmpi.dll msmpi.dll
  #     gendef msmpi.dll
  #     dlltool -d msmpi.def -l libmsmpi.a -D msmpi.dll
  #     del msmpi.def
  #     copy libmsmpi.a %MSMPI_LIB32%/libmsmpi.a

  if( WIN32 AND DEFINED CMAKE_Fortran_COMPILER AND
      TARGET MPI::MPI_Fortran )

    # only do this if we are in a CMakeAddFortranSubdirectory directive when
    # the main Generator is Visual Studio and the Fortran subdirectory uses
    # gfortran with Makefiles.

    # MS-MPI has an architecture specific include directory that FindMPI.cmake
    # doesn't seem to pickup correctly.  Add it here.
    get_target_property(mpilibdir MPI::MPI_Fortran
      INTERFACE_LINK_LIBRARIES)
    get_target_property(mpiincdir MPI::MPI_Fortran
      INTERFACE_INCLUDE_DIRECTORIES)
    foreach( arch x86 x64 )
      string(FIND "${mpilibdir}" "lib/${arch}" found_lib_arch)
      string(FIND "${mpiincdir}" "include/${arch}" found_inc_arch)
      if( ${found_lib_arch} GREATER 0 AND
          ${found_inc_arch} LESS 0 )
        if( IS_DIRECTORY "${mpiincdir}/${arch}")
          list(APPEND mpiincdir "${mpiincdir}/${arch}")
        endif()
      endif()
    endforeach()

    # Reset the include directories for MPI::MPI_Fortran to pull in the extra
    # $arch locations (if any)
    set_target_properties(MPI::MPI_Fortran
      PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${mpiincdir}")

  endif()

  if( ${MPI_C_FOUND} )
    message(STATUS "Looking for MPI.......found ${MPIEXEC_EXECUTABLE}")
  else()
    message(STATUS "Looking for MPI.......not found")
  endif()

  mark_as_advanced( MPI_FLAVOR MPIEXEC_OMP_POSTFLAGS MPI_LIBRARIES )

endmacro( setupMPILibrariesWindows )

#----------------------------------------------------------------------#
# End setupMPI.cmake
#----------------------------------------------------------------------#
