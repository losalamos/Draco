#-----------------------------*-cmake-*----------------------------------------#
# file   draco_regression_macros.cmake
# brief  Helper macros for setting up a CTest/CDash regression system
# note   Copyright (C) 2016-2018 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

# Ref: http://www.cmake.org/Wiki/CMake_Testing_With_CTest
#      http://www.cmake.org/Wiki/CMake_Scripting_of_CTest

# Echo settings if 'ON'
set( drm_verbose ON )

# Call this script from regress/Draco_*.cmake

##----------------------------------------------------------------------------##
## Testing parallelism:
## Use the result of this command in the main ctest script like this:
##     find_num_procs_avail_for_running_tests( num_test_procs )
##     set(ctest_test_args ${ctest_test_args} PARALLEL_LEVEL ${num_test_procs})
## Find the number of processors that can be used for testing.
##----------------------------------------------------------------------------##
macro( find_num_procs_avail_for_running_tests )

  # Default value is 1
  unset( num_test_procs )
  unset( max_system_load )

  # If this job is running under Torque (msub script), use the environment
  # variable PBS_NP or SLURM_NPROCS
  if( DEFINED ENV{CRAYPE_DIR} )

    set( num_test_procs 1 )
    set( num_sockets 1 )
    if( EXISTS /usr/bin/lscpu )
      execute_process( COMMAND /usr/bin/lscpu
        OUTPUT_VARIABLE lscpu_out
        OUTPUT_STRIP_TRAILING_WHITESPACE )
      # break text block into a list of lines
      string( REPLACE "\n" ";"  lscpu ${lscpu_out} )
      foreach( line ${lscpu} )
        if( ${line} MATCHES "Socket")
          string( REGEX REPLACE "^.* ([0-9]+)$" "\\1" num_sockets ${line} )
        endif()
        # find number of physical cores
        if( ${line} MATCHES "per socket")
          string( REGEX REPLACE "^.* ([0-9]+)$" "\\1" num_test_procs ${line} )
        endif()
        # if( ${line} MATCHES "Model name")
        #   if( ${line} MATCHES "Xeon Phi" )
        #     #     message("is knl")
        #     # With srun, the KNL has trouble with 68 cores worth of jobs, so
        #     # limit to a lower number
        #     set( num_test_procs 8 )
        #   endif()
        # endif()
      endforeach()
      math( EXPR num_test_procs "${num_test_procs} * ${num_sockets}" )
      unset( num_sockets )
    endif()

  elseif( NOT "$ENV{PBS_NP}x" STREQUAL "x" )
    set( num_test_procs $ENV{PBS_NP} )
  elseif( NOT "$ENV{SLURM_NPROCS}x" STREQUAL "x")
    set( num_test_procs $ENV{SLURM_NPROCS} )
  else()
    # If this is not a known batch system, then attempt to set values according
    # to machine properties
    include(ProcessorCount)
    ProcessorCount(num_test_procs)
  endif()

  # Limit the number of concurrent tests running by looking at the system load:
  math( EXPR max_system_load "${num_test_procs} + 1" )

endmacro()

# ------------------------------------------------------------
# Defaults (override with optional arguments)
# ------------------------------------------------------------
macro( set_defaults )

  message("\n----------------------------------------
Setting defaults
----------------------------------------\n")

  # Prerequisites:
  #
  # This setup assumes that the project work_dir will contain 3 subdirectories:
  # source, build and target.  See how CMAKE_SOURCE_DIRECTORY,
  # CMAKE_BUILD_DIRECTORY AND CMAKE_INSTALL_PREFIX are set just below.
  if( NOT work_dir )
    if( EXISTS "$ENV{work_dir}" )
      set( work_dir $ENV{work_dir} )
    else( EXISTS "$ENV{work_dir}" )
      if( EXISTS "$ENV{PWD}" )
        set( work_dir $ENV{PWD} )
        message( "
Warning:  work_dir has been set to pwd.
Set work_dir in your environment if you want to use a different
location.  You must also ensure that work_dir exists!
" )
      else( EXISTS "$ENV{PWD}" )
        message( FATAL_ERROR "
You must set the variable \"work_dir\" either as an environment
variable or on the ctest command line.  The value of the variable
should be an existing directory path where you want the regression
build to take place. Try:
unix% export work_dir=/full/path/to/work_dir
or
win32$ set work_dir=c:/full/path/to/work_dir
")
      endif( EXISTS "$ENV{PWD}" )
    endif( EXISTS "$ENV{work_dir}" )
  endif( NOT work_dir )
  file( TO_CMAKE_PATH ${work_dir} work_dir )

  # Set the sitename, but strip any domain information. If we are on an HPC
  # machine, attempt to associate the backend name with the same string that is
  # used to identify the front end. If we are on a LANL HPC machine, attempt to
  # use the sys_name tool for this purpose.
  if( EXISTS /usr/projects/hpcsoft/utilities/bin/sys_name )
    execute_process( COMMAND /usr/projects/hpcsoft/utilities/bin/sys_name
      OUTPUT_VARIABLE sitename
      OUTPUT_STRIP_TRAILING_WHITESPACE )
  else()
    site_name( sitename )
    string( REGEX REPLACE "([A-z0-9]+).*" "\\1" sitename ${sitename} )
  endif()

  if( ${sitename} MATCHES "cn[0-9]+" OR ${sitename} MATCHES "darwin-fe")
    set( sitename "Darwin" )
  endif()
  message( "sitename = ${sitename}")
  set( CTEST_SITE ${sitename} )
  set( CTEST_SOURCE_DIRECTORY "${work_dir}/source" )
  set( CTEST_BINARY_DIRECTORY "${work_dir}/build"  )
  set( CMAKE_INSTALL_PREFIX   "${work_dir}/target" )

  # Default is "Experimental." Special builds are "Nightly" or "Continuous"
  set( CTEST_MODEL "Experimental" )

  # Default is "Release."
  # Special types are "Debug," "RelWithDebInfo" or "MinSizeRel"
  set( CTEST_BUILD_CONFIGURATION "Release" )

  if( WIN32 )
    # add option for "NMake Makefiles JOM"?
    # set( CTEST_CMAKE_GENERATOR "NMake Makefiles" )
    set( CTEST_CMAKE_GENERATOR "Visual Studio 15 2017 Win64" )
  else()
    set( CTEST_CMAKE_GENERATOR "Unix Makefiles" )
  endif()

  set( CTEST_USE_LAUNCHERS 0 )
  set( ENABLE_C_CODECOVERAGE OFF )

  # Dashboard setup (in place of CTestConfig.cmake)
  if( NOT CTEST_PROJECT_NAME )
     set( CTEST_PROJECT_NAME "UnknownProject")
  endif()
  # ALL CRON JOBS MUST START AFTER THIS TIME + 1 HOUR (FOR DST).
  # This should be set in each projects CTestConfig.cmake file.
  #set( CTEST_NIGHTLY_START_TIME "00:00:01 MST")

  set( CTEST_DROP_METHOD "https")
  # if( "${CTEST_SITE}" MATCHES "Darwin" )
  #   set( CTEST_DROP_METHOD "http")
  # endif()
  set( CTEST_DROP_SITE "rtt.lanl.gov")
  set( CTEST_DROP_LOCATION
     "/cdash/submit.php?project=${CTEST_PROJECT_NAME}" )
  set( CTEST_DROP_SITE_CDASH TRUE )
  set( CTEST_CURL_OPTIONS CURLOPT_SSL_VERIFYPEER_OFF CURLOPT_SSL_VERIFYHOST_OFF )

  if( EXISTS "$ENV{VENDOR_DIR}" )
    file( TO_CMAKE_PATH $ENV{VENDOR_DIR} VENDOR_DIR )
  endif()
  find_path( VENDOR_DIR
    ChangeLog
    PATHS
      /scratch/vendors
      /usr/projects/draco/vendors
      c:/vendors
      )
   if( NOT "$ENV{AUTODOCDIR}x" STREQUAL "x" )
     set( AUTODOCDIR "$ENV{AUTODOCDIR}" )
     get_filename_component( AUTODOCDIR "${AUTODOCDIR}" ABSOLUTE )
   endif()

   set( VERBOSE ON )
   set( CTEST_OUTPUT_ON_FAILURE ON )

   # The default timeout is 10 min, change this to 30 min.
   set( CTEST_TEST_TIMEOUT "1800" ) # seconds

   # Parallelism for the regression:
   # - Assume that the compile phase done on the same platform/partition as the
   #   configure phase.  Query the system using ProcessorCount to find out how
   #   many cores we can use.
   include(ProcessorCount)
   ProcessorCount(num_compile_procs)
   if( NOT WIN32 )
       if(NOT num_compile_procs EQUAL 0)
         set(CTEST_BUILD_FLAGS "-j${num_compile_procs} -l${num_compile_procs}")
         if( "${sitename}" STREQUAL "Trinity" OR "${sitename}" STREQUAL "Trinitite")
           # We compile on the front end for this machine. Since we don't know
           # the actual load a priori, we use the -l option to limit the total
           # load on the machine.  For CT, my experience shows that 'make -l N'
           # actually produces a machine load ~ 1.5*N, so we will specify the
           # max load to be half of the total number of procs.
           math(EXPR half_num_compile_procs "${num_compile_procs} / 2" )
           set(CTEST_BUILD_FLAGS "-j ${half_num_compile_procs} -l ${num_compile_procs}")
         endif()
       endif()
   else()
     if(NOT num_compile_procs EQUAL 0)
       # Parallel builds for 'msbuild'
       # Ref: https://docs.microsoft.com/en-us/visualstudio/msbuild/msbuild-command-line-reference?view=vs-2015
       set(CTEST_BUILD_FLAGS "-m:${num_compile_procs}")
     endif()
   endif()

   # Testing parallelism
   # - Call the macro
   #   find_num_procs_avail_for_running_tests(num_test_procs) from the each
   #   projects ctest script (Draco_Linux64.ctest).

   # Echo settings
   if( ${drm_verbose} )
     message("
work_dir   = ${work_dir}

CTEST_PROJECT_NAME     = ${CTEST_PROJECT_NAME}
CTEST_SCRIPT_DIRECTORY = ${CTEST_SCRIPT_DIRECTORY}
CTEST_SCRIPT_NAME      = ${CTEST_SCRIPT_NAME}

CTEST_SITE             = ${CTEST_SITE}
CTEST_SOURCE_DIRECTORY = ${CTEST_SOURCE_DIRECTORY}
CTEST_BINARY_DIRECTORY = ${CTEST_BINARY_DIRECTORY}
CMAKE_INSTALL_PREFIX   = ${CMAKE_INSTALL_PREFIX}
CTEST_MODEL            = ${CTEST_MODEL}
CTEST_BUILD_CONFIGURATION = ${CTEST_BUILD_CONFIGURATION}
CTEST_CMAKE_GENERATOR  = ${CTEST_CMAKE_GENERATOR}
CTEST_NIGHTLY_START_TIME  = ${CTEST_NIGHTLY_START_TIME}
CTEST_DROP_METHOD         = ${CTEST_DROP_METHOD}
CTEST_DROP_SITE           = ${CTEST_DROP_SITE}
CTEST_DROP_LOCATION       = ${CTEST_DROP_LOCATION}
CTEST_DROP_SITE_CDASH     = ${CTEST_DROP_SITE_CDASH}
CTEST_CURL_OPTIONS        = ${CTEST_CURL_OPTIONS}
CTEST_BUILD_FLAGS         = ${CTEST_BUILD_FLAGS}
num_compile_procs         = ${num_compile_procs}
VENDOR_DIR                = ${VENDOR_DIR}
")
  endif()

endmacro( set_defaults )

# ------------------------------------------------------------
# Parse Arguments
# ------------------------------------------------------------
macro( parse_args )

  message("\n----------------------------------------
Parsing arguments
----------------------------------------\n")

  # Default is "Experimental." Special builds are "Nightly" or "Continuous"
  if( ${CTEST_SCRIPT_ARG} MATCHES Nightly )
    set( CTEST_MODEL "Nightly" )
  elseif( ${CTEST_SCRIPT_ARG} MATCHES Continuous  )
    set( CTEST_MODEL "Continuous" )
  endif()

  # Default is "Release."
  # Special types are "Debug," "RelWithDebInfo" or "MinSizeRel"
  if( ${CTEST_SCRIPT_ARG} MATCHES Debug )
    set( CTEST_BUILD_CONFIGURATION "Debug" )
  elseif( ${CTEST_SCRIPT_ARG} MATCHES RelWithDebInfo )
    set( CTEST_BUILD_CONFIGURATION "RelWithDebInfo" )
  elseif( ${CTEST_SCRIPT_ARG} MATCHES MinSizeRel )
    set( CTEST_BUILD_CONFIGURATION "MinSizeRel" )
  endif()

  # Post options: SubmitOnly or NoSubmit
  set( CTEST_CONFIGURE OFF )
  set( CTEST_BUILD     OFF )
  set( CTEST_TEST      OFF )
  set( CTEST_SUBMIT    OFF )
  set( CTEST_AUTODOC   OFF )
  if( ${CTEST_SCRIPT_ARG} MATCHES Configure )
    set( CTEST_CONFIGURE "ON" )
  endif()
  if( ${CTEST_SCRIPT_ARG} MATCHES Build )
    set( CTEST_BUILD "ON" )
  endif()
  if( ${CTEST_SCRIPT_ARG} MATCHES Test )
    set( CTEST_TEST "ON" )
  endif()
  if( ${CTEST_SCRIPT_ARG} MATCHES Submit )
    set( CTEST_SUBMIT "ON" )
  endif()
  if( ${CTEST_SCRIPT_ARG} MATCHES Autodoc )
    set( CTEST_AUTODOC "ON" )
  endif()

  # default compiler name based on platform
  unset(compiler_version)
  if( DEFINED ENV{LCOMPILERVER} )
    set( compiler_version $ENV{LCOMPILERVER} )
  elseif( "$ENV{CXX}" MATCHES ".*[-]([0-9]+[.][0-9]+[.-][0-9]+).*" )
    string( REGEX REPLACE ".*[-]([0-9]+[.][0-9]+[.-][0-9]+).*" "\\1"
      compiler_version "$ENV{CXX}" )
  endif()

  # refine compiler short name.
  set(USE_CUDA OFF)
  if( "$ENV{CXX}" MATCHES "pgCC" OR "$ENV{CXX}" MATCHES "pgc[+][+]" )
    set( compiler_short_name "pgi" )
  elseif("$ENV{CXX}" MATCHES "clang" )
    set( compiler_short_name "clang" )
  elseif("$ENV{CXX}" MATCHES "icpc" )
    set( compiler_short_name "intel" )
  elseif( DEFINED ENV{PE_ENV} )
    # Trinity define this variable to define the flavor of the PrgEnv module
    # currently loaded.
    if( $ENV{PE_ENV} MATCHES "INTEL")
      set( compiler_short_name "intel" )
    elseif( $ENV{PE_ENV} MATCHES "PGI")
      set( compiler_short_name "pgi" )
    elseif( $ENV{PE_ENV} MATCHES "CRAY")
      set( compiler_short_name "cray" )
    elseif( $ENV{PE_ENV} MATCHES "GNU")
      set( compiler_short_name "gcc" )
    endif()
    if( NOT DEFINED compiler_version )
      if( ${work_dir} MATCHES ".*[-]([0-9]+[.][0-9]+[.-][0-9]+).*" )
        string( REGEX REPLACE ".*[-]([0-9]+[.][0-9]+[.-][0-9]+).*" "\\1"
          compiler_version ${work_dir} )
      endif()
    endif()
  elseif( "$ENV{CC}" MATCHES "gcc" )
    set(compiler_short_name "gcc" )
    # /usr/bin/gcc
    # /ccs/codes/radtran/vendors/bullseyecoverage-8.9.75/bin/gcc
    execute_process( COMMAND $ENV{CC} --version
      OUTPUT_VARIABLE cxx_version
      OUTPUT_STRIP_TRAILING_WHITESPACE )
    string( REGEX REPLACE "[^0-9]*([0-9]+).([0-9]+).([0-9]+).*" "\\1.\\2.\\3"
      compiler_version ${cxx_version} )
  elseif( CTEST_CMAKE_GENERATOR MATCHES "Visual Studio" )
    set( compiler_short_name "cl" )
    execute_process( COMMAND cl
      ERROR_VARIABLE cxx_version
      OUTPUT_STRIP_TRAILING_WHITESPACE
      OUTPUT_QUIET
      )
    string( REGEX REPLACE "[^0-9]*([0-9]+).([0-9]+).([0-9]+).*" "\\1.\\2.\\3"
      compiler_version "${cxx_version}" )
  else()
    set( compiler_short_name "unknown" )
    set( compiler_version "unknown" )
  endif()

  # if known, append the compiler version
  set( compiler_short_name "${compiler_short_name}-${compiler_version}" )

  # append the compiler_short_name with the extra_params string (if any) and set
  # some variables based on extra_param's value.
  if( DEFINED ENV{extra_params_sort_safe} )

    if( NOT "$ENV{extra_params_sort_safe}empty" STREQUAL "empty" )
      set( compiler_short_name
        "${compiler_short_name}-$ENV{extra_params_sort_safe}" )
    endif()

    if( $ENV{extra_params_sort_safe} MATCHES "cuda" )
      set(USE_CUDA ON)
    endif()
    if( $ENV{extra_params_sort_safe} MATCHES "fulldiagnostics" )
      set( FULLDIAGNOSTICS "DRACO_DIAGNOSTICS:STRING=7")
      # Note 'DRACO_TIMING:STRING=2' will break milagro tests (python cannot
      # parse output).
    endif()
    if( $ENV{extra_params_sort_safe} MATCHES "scalar" )
      set( DRACO_C4 "DRACO_C4:STRING=SCALAR" )
    elseif( $ENV{extra_params_sort_safe} MATCHES "static" )
      set( DRACO_LIBRARY_TYPE "DRACO_LIBRARY_TYPE:STRING=STATIC" )
    endif()
    if( $ENV{extra_params_sort_safe} MATCHES "vtest" )
      list( APPEND CUSTOM_VARS "RUN_VERIFICATION_TESTS:BOOL=ON" )
    endif()
    if( $ENV{extra_params_sort_safe} MATCHES "perfbench" )
      list( APPEND CUSTOM_VARS "ENABLE_PERFBENCH:BOOL=ON" )
    endif()
  endif()

  # Set the build name: (<platform>-<compiler>-<configuration>)
  if( WIN32 )
    set( CTEST_BUILD_NAME "${compiler_short_name}-${CTEST_BUILD_CONFIGURATION}" )
    # if( "$ENV{dirext}" MATCHES "x64" )
    # endif()
  elseif( APPLE ) # OS/X
    set( CTEST_BUILD_NAME "${compiler_short_name}_${CTEST_BUILD_CONFIGURATION}" )
  else() # Unix
    set( CTEST_BUILD_NAME "${compiler_short_name}_${CTEST_BUILD_CONFIGURATION}" )
    if( "$ENV{featurebranch}notset" STREQUAL "notset" )
      message(FATAL_ERROR "Checkout from github requested, but ENV{featurebranch} is not set.")
    endif()
    set( CTEST_BUILD_NAME "${compiler_short_name}_${CTEST_BUILD_CONFIGURATION}-$ENV{featurebranch}" )
  endif()

  # I think this is only used by msbuild targets (e.g. 'ctest -j 4 -C Debug'),
  # but I want to define it universally.
  set( CTEST_CONFIGURATION_TYPE "${CTEST_BUILD_CONFIGURATION}" )

  # Default is no Coverage Analysis
  if( ${CTEST_SCRIPT_ARG} MATCHES Coverage )
    if( ${CTEST_BUILD_CONFIGURATION} MATCHES Release OR
        ${CTEST_BUILD_CONFIGURATION} MATCHES MinSizeRel )
      message( FATAL_ERROR "Cannot run coverage for \"Release\" mode builds." )
    endif()
    set( ENABLE_C_CODECOVERAGE ON )
    set( CTEST_BUILD_NAME "${CTEST_BUILD_NAME}_Cov" )
  endif()

  # DynamicAnalysis (valgrind)
  if( ${CTEST_SCRIPT_ARG} MATCHES DynamicAnalysis )
    set( ENABLE_DYNAMICANALYSIS ON )
    set( CTEST_BUILD_NAME "${CTEST_BUILD_NAME}_DA" )
  endif()

  # Bounds Checking
  if( ${CTEST_SCRIPT_ARG} MATCHES bounds_checking )
    if( "${compiler_short_name}" MATCHES "gcc-4.8.5" )
      set( BOUNDS_CHECKING "GCC_ENABLE_GLIBCXX_DEBUG:BOOL=ON" )
    else()
      message(FATAL_ERROR "I don't know how to turn on bounds checking for compiler = ${compiler_short_name}" )
    endif()
  endif()

  if( DEFINED ENV{buildname_append} )
    set( CTEST_BUILD_NAME "${CTEST_BUILD_NAME}$ENV{buildname_append}" )
  endif()

# Convert CUSTOM_VARS from a list into a multi-line string
string(REGEX REPLACE ";" "\n" CUSTOM_VARS "${CUSTOM_VARS}")

  if( ${drm_verbose} )
    message("
CTEST_MODEL                 = ${CTEST_MODEL}
CTEST_BUILD_CONFIGURATION   = ${CTEST_BUILD_CONFIGURATION}
CTEST_CONFIGURATION_TYPE    = ${CTEST_CONFIGURATION_TYPE}
compiler_short_name         = ${compiler_short_name}
compiler_version            = ${compiler_version}
CTEST_BUILD_NAME            = ${CTEST_BUILD_NAME}
ENABLE_C_CODECOVERAGE       = ${ENABLE_C_CODECOVERAGE}
ENABLE_DYNAMICANALYSIS      = ${ENABLE_DYNAMICANALYSIS}
CTEST_USE_LAUNCHERS         = ${CTEST_USE_LAUNCHERS}
CUSTOM_VARS                 = ${CUSTOM_VARS}
")
  endif()
endmacro( parse_args )

# ------------------------------------------------------------
# Names of developer tools
# ------------------------------------------------------------
macro( find_tools )

  message("\n----------------------------------------
Finding tools...
----------------------------------------\n")

  find_program( CTEST_CMD
    NAMES ctest
    HINTS
      "c:/Program Files (x86)/CMake 3.0/bin"
      # NO_DEFAULT_PATH
    )
  if( NOT EXISTS ${CTEST_CMD} )
    message( FATAL_ERROR "Cound not find ctest executable.(CTEST_CMD = ${CTEST_CMD})" )
  endif()

  find_program( CTEST_GIT_COMMAND
     NAMES git
     HINTS
        "C:/Program Files/Git/bin"
        # NO_DEFAULT_PATH
     )
  if( NOT EXISTS "${CTEST_GIT_COMMAND}" )
    message( FATAL_ERROR "Cound not find git executable." )
  endif()

  find_program( CTEST_CMAKE_COMMAND
    NAMES cmake
    HINTS
      "c:/Program Files (x86)/CMake 3.0/bin"
      # NO_DEFAULT_PATH
    )
  if( NOT EXISTS "${CTEST_CMAKE_COMMAND}" )
    message( FATAL_ERROR "Cound not find cmake executable." )
  endif()

  if( NOT WIN32 )
    # if MAKECOMMAND is found when using "Visual Studio" as the generator,
    # the compiler 'cl' will be found to be unable to compile a simple
    # program.
    find_program( MAKECOMMAND NAMES make )
    # No memory check program on Windows for now.
    find_program( CTEST_MEMORYCHECK_COMMAND NAMES valgrind )
  endif()

  if(ENABLE_C_CODECOVERAGE)
    find_program( COV01 NAMES cov01 )
    if( COV01 )
       get_filename_component( beyedir ${COV01} PATH )
       set( CC ${beyedir}/gcc )
       set( CXX ${beyedir}/g++ )
       set( ENV{CC} ${beyedir}/gcc )
       set( ENV{CXX} ${beyedir}/g++ )
       set( RES 1 )
       execute_process(COMMAND ${COV01} --on RESULT_VARIABLE RES )
       if( RES )
          message(FATAL_ERROR "could not run cov01 -1")
       else()
          message("BullseyeCoverage turned on")
       endif()
    else()
       message( FATAL_ERROR "Coverage requested, but bullseyecoverage's cov01 binary not in PATH." )
    endif()
  endif()

  # This breaks NMake Makefile builds because it is missing the -G"..."
  # set(CTEST_CONFIGURE_COMMAND "${CMAKE_COMMAND} \"${CTEST_SOURCE_DIRECTORY}\"")
  if( ${drm_verbose} )
    message("
CTEST_CMD           = ${CTEST_CMD}
CTEST_GIT_COMMAND   = ${CTEST_GIT_COMMAND}
CTEST_CMAKE_COMMAND = ${CTEST_CMAKE_COMMAND}
MAKECOMMAND         = ${MAKECOMMAND}
CTEST_MEMORYCHECK_COMMAND         = ${CTEST_MEMORYCHECK_COMMAND}
MEMORYCHECK_SUPPRESSIONS_FILE     = ${MEMORYCHECK_SUPPRESSIONS_FILE}
CTEST_MEMORYCHECK_COMMAND_OPTIONS = ${CTEST_MEMORYCHECK_COMMAND_OPTIONS}

")
    if(ENABLE_C_CODECOVERAGE)
       message("
CC    = ${CC}
CXX   = ${CXX}
COV01 = ${COV01}
RES   = ${RES}
")
    endif()
  endif()
endmacro( find_tools )

# ------------------------------------------------------------
# Setup SVNROOT
# ------------------------------------------------------------
macro( set_svn_command svnpath )
  if( NOT EXISTS ${CTEST_SOURCE_DIRECTORY}/CMakeLists.txt )
    if( EXISTS /ccs/codes/radtran/svn )
      set( CTEST_CHECKOUT_COMMAND
        "${CTEST_CVS_COMMAND} checkout file:///ccs/codes/radtran/svn/${svnpath} source" )
      message("CTEST_CHECKOUT_COMMAND = ${CTEST_CHECKOUT_COMMAND}")
    elseif( EXISTS /usr/projects/draco/regress/svn ) # CCS-7's Darwin
      set( CTEST_CHECKOUT_COMMAND
        "${CTEST_CVS_COMMAND} checkout file:///usr/projects/draco/regress/svn/${svnpath} source" )
    elseif( EXISTS /usr/projects/jayenne/regress/svn ) # HPC machines (ML, CT)
      set( CTEST_CHECKOUT_COMMAND
        "${CTEST_CVS_COMMAND} checkout file:///usr/projects/jayenne/regress/svn/${svnpath} source" )
    else()
      set( CTEST_CHECKOUT_COMMAND
        "${CTEST_CVS_COMMAND} checkout svn+ssh://ccscs7/ccs/codes/radtran/svn/${svnpath} source" )
    endif()
  endif()
endmacro()

macro( set_git_command gitpath )
  set( CTEST_UPDATE_TYPE "git" )
  if( NOT EXISTS ${CTEST_SOURCE_DIRECTORY}/CMakeLists.txt )
    if( ${gitpath} MATCHES "Draco" )
      set( CTEST_CHECKOUT_COMMAND
        "${CTEST_GIT_COMMAND} clone --depth 1 https://github.com/lanl/${gitpath} source" )
    else()
      # This assumes that a valid ssh-key exists in the current environment and
      # works with gitlab.lanl.gov.
      if( IS_DIRECTORY "$ENV{gitroot}" )
        set( gitroot "$ENV{gitroot}/" )
        set( CTEST_CHECKOUT_COMMAND
          "${CTEST_GIT_COMMAND} clone --depth 1 ${gitroot}${gitpath} source" )
      else()
        set( gitroot "git@gitlab.lanl.gov:" )
        set( CTEST_CHECKOUT_COMMAND
          "${CTEST_GIT_COMMAND} clone --depth 1 ${gitroot}${gitpath} source" )
      endif()
    endif()
  endif()
  # normaly, just use the 'develop' branch.  Otherwise ENV{featurebranch} will
  # be set to something like pr42.
  if( "$ENV{featurebranch}notset" STREQUAL "notset" OR
      "$ENV{featurebranch}" STREQUAL "develop" )
    set( CTEST_GIT_UPDATE_CUSTOM "${CTEST_GIT_COMMAND};pull;origin;develop")
  elseif( "$ENV{featurebranch}" MATCHES "pr[0-9]+" )
    string( REPLACE "pr" "" featurebranch "$ENV{featurebranch}" )
    # set( CTEST_GIT_UPDATE_CUSTOM "${CTEST_GIT_COMMAND};fetch;origin;pull/${featurebranch}/head:pr${featurebranch};&&;${CTEST_GIT_COMMAND} checkout pr${featurebranch}")
    # assumes that this is run when PWD = ${CTEST_SOURCE_DIRECTORY} and git
    # clone was successful.
    # set( CTEST_GIT_UPDATE_CUSTOM "${CTEST_SOURCE_DIRECTORY}/regression/fetch_co.sh;${CTEST_GIT_COMMAND};${featurebranch}")
    if ( "${githost}x" STREQUAL "x" )
      if( ${gitpath} MATCHES "Draco" )
        set( githost "github")
      else()
        set( githost "gitlab")
      endif()
    endif()
    set( CTEST_GIT_UPDATE_CUSTOM "$ENV{rscriptdir}/fetch_co.sh;${CTEST_GIT_COMMAND};${githost};${featurebranch}")
  else()
    message( FATAL_ERROR "I don't know how to checkout git feature branch named '$ENV{featurebranch}'.")
  endif()

  message( "
CTEST_UPDATE_TYPE       = ${CTEST_UPDATE_TYPE}
CTEST_CHECKOUT_COMMAND  = ${CTEST_CHECKOUT_COMMAND}
CTEST_GIT_UPDATE_CUSTOM = ${CTEST_GIT_UPDATE_CUSTOM}
")

endmacro()

# ------------------------------------------------------------
# Setup for Code Coverage and LOC metrics
# ------------------------------------------------------------
macro( setup_for_code_coverage )
   if( "${sitename}" MATCHES "ccscs[1-9]" AND "$ENV{CXX}" MATCHES "g[+][+]" )
      if( ${CTEST_BUILD_CONFIGURATION} MATCHES Debug )
         if(ENABLE_C_CODECOVERAGE)

            # Code coverage setup
            message("Generating ${CTEST_BINARY_DIRECTORY}/covclass_cmake.cfg")
            configure_file(
               ${CTEST_SCRIPT_DIRECTORY}/covclass_cmake.cfg
               ${CTEST_BINARY_DIRECTORY}/covclass_cmake.cfg
               @ONLY )
            set( ENV{COVDIRCFG}   ${CTEST_BINARY_DIRECTORY}/covclass_cmake.cfg )
            set( ENV{COVFNCFG}    ${CTEST_BINARY_DIRECTORY}/covclass_cmake.cfg )
            set( ENV{COVCLASSCFG} ${CTEST_BINARY_DIRECTORY}/covclass_cmake.cfg )
            set( ENV{COVSRCCFG}   ${CTEST_BINARY_DIRECTORY}/covclass_cmake.cfg )
            set( ENV{COVFILE}     ${CTEST_BINARY_DIRECTORY}/CMake.cov )
            execute_process(COMMAND "${COV01}" --on --verbose RESULT_VARIABLE RES)

            # Process and save lines of code
            message( "Generating lines of code statistics...
/scratch/vendors/bin/cloc
               --exclude-dir=heterogeneous,chimpy
               --exclude-list-file=$ENV{rscriptdir}/cloc-exclude.cfg
               --exclude-lang=Text,Postscript
               --categorize=cloc-categorize.log
               --counted=cloc-counted.log
               --ignored=cloc-ignored.log
               --progress-rate=0
               --report-file=lines-of-code.log
               --force-lang-def=$ENV{rscriptdir}/cloc-lang.defs
               ${CTEST_SOURCE_DIRECTORY}
            ")
            execute_process(
               COMMAND /scratch/vendors/bin/cloc
               --exclude-dir=heterogeneous,chimpy
               --exclude-list-file=$ENV{rscriptdir}/cloc-exclude.cfg
               --exclude-lang=Text,Postscript
               --categorize=cloc-categorize.log
               --counted=cloc-counted.log
               --ignored=cloc-ignored.log
               --progress-rate=0
               --report-file=lines-of-code.log
               --force-lang-def=$ENV{rscriptdir}/cloc-lang.defs
               ${CTEST_SOURCE_DIRECTORY}
               #  --3
               #  --diff
               WORKING_DIRECTORY ${CTEST_BINARY_DIRECTORY}
               )
            message( "Lines of code data at ${CTEST_BINARY_DIRECTORY}/lines-of-code.log")
            message( "Generating lines of code statistics (omitting test directories)
/scratch/vendors/bin/cloc
               --exclude-dir=heterogeneous,chimpy,test
               --exclude-list-file=/scratch/regress/draco/regression/cloc-exclude.cfg
               --exclude-lang=Text,Postscript
               --categorize=cloc-categorize-notest.log
               --counted=cloc-counted-notest.log
               --ignored=cloc-ignored-notest.log
               --progress-rate=0
               --report-file=lines-of-code-notest.log
               --force-lang-def=/scratch/regress/draco/regression/cloc-lang.defs
               ${CTEST_SOURCE_DIRECTORY}
            ")
            execute_process(
               COMMAND /scratch/vendors/bin/cloc
               --exclude-dir=heterogeneous,chimpy,test
               --exclude-list-file=/scratch/regress/draco/regression/cloc-exclude.cfg
               --exclude-lang=Text,Postscript
               --categorize=cloc-categorize.log
               --counted=cloc-counted.log
               --ignored=cloc-ignored.log
               --progress-rate=0
               --report-file=lines-of-code-notest.log
               --force-lang-def=/scratch/regress/draco/regression/cloc-lang.defs
               ${CTEST_SOURCE_DIRECTORY}
               #  --3
               #  --diff
               WORKING_DIRECTORY ${CTEST_BINARY_DIRECTORY}
               )
            message( "Lines of code data at ${CTEST_BINARY_DIRECTORY}/lines-of-code.log")

            # set( CTEST_NOTES_FILES file1 file2 )

            # create a cdash link to redmine
            # http://www.kitware.com/source/home/post/59
            # https://github.com/Slicer/Slicer/blob/57f14d0d233ee103e365161cfc0b3962df0bc203/CMake/MIDASCTestUploadURL.cmake#L69
            # In CDash, ensure that
            # Settings->Project->Miscellaneous->File upload quota > 0.
            #file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/redmine.url"
            #   "http://rtt.lanl.gov/redmine" )
            ctest_upload( FILES
            #  "${CMAKE_CURRENT_BINARY_DIR}/redmine.url"
              "${CTEST_BINARY_DIRECTORY}/lines-of-code.log"
              "${CTEST_BINARY_DIRECTORY}/lines-of-code-notest.log" )

         endif()
      endif()
   endif()
endmacro( setup_for_code_coverage )

# ------------------------------------------------------------
# Process Code Coverage or Dynamic Memory Analysis
#
# If BUILD_CONFIG = Debug do dynamic analysis (valgrind)
# If BUILD_CONFIG = Coverage to code coverage (bullseye)
#
# dyanmic analysis excludes tests with label "nomemcheck"
# ------------------------------------------------------------
macro(process_cc_or_da)
  if(ENABLE_C_CODECOVERAGE)
    message( "ctest_coverage( BUILD \"${CTEST_BINARY_DIRECTORY}\" )")
    ctest_coverage( BUILD "${CTEST_BINARY_DIRECTORY}" )
    message( "Generating code coverage log file: ${CTEST_BINARY_DIRECTORY}/covdir.log

cd ${CTEST_BINARY_DIRECTORY}
covdir -o ${CTEST_BINARY_DIRECTORY}/covdir.log")
    execute_process(
      COMMAND
      covdir -o ${CTEST_BINARY_DIRECTORY}/covdir.log
      WORKING_DIRECTORY ${CTEST_BINARY_DIRECTORY}
      )
    list( APPEND CTEST_NOTES_FILES "${CTEST_BINARY_DIRECTORY}/covdir.log")
    execute_process(COMMAND "${COV01}" --off RESULT_VARIABLE RES)
  endif()

  if( ${ENABLE_DYNAMICANALYSIS} )
    set( CTEST_TEST_TIMEOUT 7200 )
    if( EXISTS ${CTEST_SOURCE_DIRECTORY}/regression/valgrind_suppress.txt )
      set( MEMORYCHECK_SUPPRESSIONS_FILE
        ${CTEST_SOURCE_DIRECTORY}/regression/valgrind_suppress.txt )
    elseif( EXISTS ${CTEST_SOURCE_DIRECTORY}/scripts/valgrind_suppress.txt )
      set( MEMORYCHECK_SUPPRESSIONS_FILE
        ${CTEST_SOURCE_DIRECTORY}/scripts/valgrind_suppress.txt )
    elseif( EXISTS ${CTEST_SCRIPT_DIRECTORY}/valgrind_suppress.txt )
      set( MEMORYCHECK_SUPPRESSIONS_FILE
        ${CTEST_SCRIPT_DIRECTORY}/valgrind_suppress.txt )
    endif()
    if( EXISTS ${MEMORYCHECK_SUPPRESSIONS_FILE} )
      set( valgrind_suppress_option
        "--suppressions=${MEMORYCHECK_SUPPRESSIONS_FILE}")
      set( CTEST_MEMORYCHECK_COMMAND_OPTIONS
        "-q --tool=memcheck --leak-check=full --trace-children=yes --gen-suppressions=all ${valgrind_suppress_option}" )
    endif()
    message(" CTEST_MEMORYCHECK_COMMAND_OPTIONS = ${CTEST_MEMORYCHECK_COMMAND_OPTIONS}")

    set( ctest_memcheck_options "SCHEDULE_RANDOM ON" )
    string( APPEND ctest_memcheck_options " PARALLEL_LEVEL ${num_test_procs}" )
    string( APPEND ctest_memcheck_options " EXCLUDE_LABEL nomemcheck")

    # if we are running on a machine that openly shares resources, use the
    # TEST_LOAD feature to limit the number of cores used while testing. For
    # machines that run schedulers, the whole allocation is available so there
    # is no need to limit the load.
    if( "${CTEST_SITE}" MATCHES "ccscs" )
      string( APPEND ctest_memcheck_options " TEST_LOAD ${max_system_load}" )
    endif()
    message( "ctest_memcheck( ${ctest_memcheck_options} )" )
    separate_arguments( ctest_memcheck_options )
    ctest_memcheck( ${ctest_memcheck_options} )
  endif()
endmacro(process_cc_or_da)

# ------------------------------------------------------------
# Special default settings for a couple of platforms
#
# Sets DRACO_DIR
# ------------------------------------------------------------
macro(set_pkg_work_dir this_pkg dep_pkg)

  string( TOUPPER ${dep_pkg} dep_pkg_caps )
  # Assume that draco_work_dir is parallel to our current location, but only
  # replace the directory name preceeding the dashboard name.
  file( TO_CMAKE_PATH "$ENV{work_dir}" work_dir )
  file( TO_CMAKE_PATH "$ENV{DRACO_DIR}" DRACO_DIR )
  string( REGEX REPLACE "${this_pkg}[/\\](Nightly|Experimental|Continuous)"
    "${dep_pkg}/\\1" ${dep_pkg}_work_dir ${work_dir} )

  # If this is a special build, link to the normal Debug/Release Draco files:
  #
  # J/C directory       Draco directory
  # --------------      -----------------
  # *-nr-*              *-*
  # *-vtest-*           *-*
  # *-perfbench-*       *-*
  # *-knl-perfbench-*   *-knl-*

  if( "${dep_pkg}" MATCHES "draco" )
    # not any ${extraparam} since many map to draco builds that have the same
    # ${extraparam}.  For example: *-newtools-*.
    foreach( extraparam nr perfbench vtest )
      string( REGEX REPLACE "[-_]${extraparam}[-_]" "-" ${dep_pkg}_work_dir
        ${${dep_pkg}_work_dir} )
    endforeach()
    if( "${this_pkg}" MATCHES "jayenne" OR "${this_pkg}" MATCHES "capsaicin" OR "${this_pkg}" MATCHES "core")
      # If this is jayenne, we might be building a pull request. Replace the PR
      # number in the path with '-develop' before looking for draco.
      string( REGEX REPLACE "(Nightly|Experimental|Continuous)_(.*)(-pr[0-9]+)/"
        "\\1_\\2-develop/" ${dep_pkg}_work_dir ${${dep_pkg}_work_dir} )
    endif()
  endif()

  find_file( ${dep_pkg}_target_dir
    NAMES README.${dep_pkg} README.md
    HINTS
      # if DRACO_DIR is defined, use it.
      ${DRACO_DIR}
      # Try a path parallel to the work_dir
      ${${dep_pkg}_work_dir}/target
    NO_DEFAULT_PATH
  )
  if( NOT EXISTS ${${dep_pkg}_target_dir} )
    message( FATAL_ERROR
      "Could not locate the ${dep_pkg} installation directory.
      ${dep_pkg}_target_dir = ${${dep_pkg}_target_dir}
      ${dep_pkg}_work_dir   = ${${dep_pkg}_work_dir}")
  endif()

  get_filename_component( ${dep_pkg_caps}_DIR ${${dep_pkg}_target_dir} PATH )
  unset( dep_pkg_caps )

endmacro(set_pkg_work_dir)

#------------------------------------------------------------------------------#
# End draco_regression_macros.cmake
#------------------------------------------------------------------------------#
