#-----------------------------*-cmake-*----------------------------------------#
# file   draco_regression_macros.cmake
# brief  Helper macros for setting up a CTest/CDash regression system
# note   Copyright (C) 2010-2012 Los Alamos National Security
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

# Ref: http://www.cmake.org/Wiki/CMake_Testing_With_CTest
#      http://www.cmake.org/Wiki/CMake_Scripting_of_CTest


# Echo settings if 'ON'
set( drm_verbose OFF )

# Call this script from regress/Draco_*.cmake

# ------------------------------------------------------------
# Defaults (override with optional arguments)
# ------------------------------------------------------------
macro( set_defaults )

  # Prerequisits:
  # 
  # This setup assumes that the project work_dir will contain 3
  # subdirectories: source, build and target.  See how 
  # CMAKE_SOURCE_DIRECTORY, CMAKE_BUILD_DIRECTORY AND
  # CMAKE_INSTALL_PREFIX are set just below.
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

  # Set the sitename, but strip any domain information
  site_name( sitename )
  string( REGEX REPLACE "([A-z0-9]+).*" "\\1" sitename ${sitename} )
  # message( "sitename = ${sitename}")
  if( ${sitename} MATCHES "ct" )
     set( sitename "Cielito" )
  elseif( ${sitename} MATCHES "ml[0-9]+" OR ${sitename} MATCHES "ml-fey")
     set( sitename "Moonlight" )
  endif()
#  message( "sitename = ${sitename}")
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
    set( CTEST_CMAKE_GENERATOR "NMake Makefiles" )
    # set( CTEST_CMAKE_GENERATOR "Visual Studio 11" )
  else()
    set( CTEST_CMAKE_GENERATOR "Unix Makefiles" )
  endif()      

  set( CTEST_USE_LAUNCHERS 0 )
  set( ENABLE_C_CODECOVERAGE OFF )
  set( ENABLE_Fortran_CODECOVERAGE OFF )

  # Dashboard setup (in place of CTestConfig.cmake)
  if( NOT CTEST_PROJECT_NAME )
     set( CTEST_PROJECT_NAME "UnknownProject")
  endif()
  # ALL CRON JOBS MUST START AFTER THIS TIME + 1 HOUR (FOR DST).
  # This should be set in each projects CTestConfig.cmake file.
  #set( CTEST_NIGHTLY_START_TIME "00:00:01 MST")
  
  set( CTEST_DROP_METHOD "https")
  set( CTEST_DROP_SITE "rtt.lanl.gov")
  set( CTEST_DROP_LOCATION 
     "/cdash/submit.php?project=${CTEST_PROJECT_NAME}" )
  set( CTEST_DROP_SITE_CDASH TRUE )
  set( CTEST_CURL_OPTIONS CURLOPT_SSL_VERIFYPEER_OFF CURLOPT_SSL_VERIFYHOST_OFF )

  set( MPIEXEC_MAX_NUMPROCS 4 CACHE STRING  "Number of cores on the local machine." )

  if( EXISTS "$ENV{VENDOR_DIR}" )
    file( TO_CMAKE_PATH $ENV{VENDOR_DIR} VENDOR_DIR )
  endif()
  find_path( VENDOR_DIR
    ChangeLog
    PATHS
      /ccs/codes/radtran/vendors/Linux64
      /usr/projects/draco/vendors/${CMAKE_SYSTEM_PROCESSOR}-${CMAKE_SYSTEM_NAME}
      /usr/projects/draco/vendors
      c:/vendors/${CMAKE_SYSTEM_PROCESSOR}-${CMAKE_SYSTEM_NAME}
      c:/vendors
      )
   set( AUTODOCDIR "${VENDOR_DIR}/../autodoc" )
   get_filename_component( AUTODOCDIR "${AUTODOCDIR}" ABSOLUTE )

   set( VERBOSE ON )
   set( CTEST_OUTPUT_ON_FAILURE ON )

   # The default timeout is 10 min, change this to 30 min.
   set( CTEST_TEST_TIMEOUT "1800" ) # seconds

  # Echo settings
  if( ${drm_verbose} )  
     message("
ARGV     = ${ARGV}

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
VENDOR_DIR                = ${VENDOR_DIR}
")
  endif()

endmacro( set_defaults )

# ------------------------------------------------------------
# Parse Arguments
# ------------------------------------------------------------
macro( parse_args )

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
  endif( ${CTEST_SCRIPT_ARG} MATCHES Debug )
  
  # Post options: SubmitOnly or NoSubmit
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

  # default compiler name based on platform
  if( WIN32 )
    set( compiler_short_name "cl" )
  else()
    set( compiler_short_name "gcc" )
  endif()
  
  # refine compiler short name.
  if( $ENV{CXX} MATCHES "pgCC" )
     set( compiler_short_name "pgi" )
  elseif($ENV{CXX} MATCHES "icpc" )
     if( ${work_dir} MATCHES ".*[-]([0-9]+[.][0-9]+[.-][0-9]+).*" )
        string( REGEX REPLACE ".*[-]([0-9]+[.][0-9]+[.-][0-9]+).*" "\\1"
           compiler_version ${work_dir} )
     elseif( ${work_dir} MATCHES ".*cuda.*" )
        set( compiler_version "cuda" )
     elseif( ${work_dir} MATCHES ".*fulldiagnostics.*" )
        set( compiler_version "fulldiagnostics" )
        set( FULLDIAGNOSTICS "DRACO_DIAGNOSTICS:STRING=7
DRACO_TIMING:STRING=2")
     endif()
     if( "${compiler_version}x" STREQUAL "x" )
        set( compiler_short_name "intel" )
     else()
        set( compiler_short_name "intel-${compiler_version}" )
     endif()
  elseif($ENV{CXX} MATCHES "xt-asyncpe" )
     # Ceilo (catamount) uses a wrapper script
     # /opt/cray/xt-asyncpe/5.06/bin/CC that masks the actual
     # compiler.  Use the following command to determine the actual
     # compiler flavor before setting compiler flags (end of this
     # macro).
     execute_process(
        COMMAND $ENV{CXX} --version
        OUTPUT_VARIABLE my_cxx_compiler
        ERROR_QUIET )
     if( ${my_cxx_compiler} MATCHES "icpc")
        set( compiler_short_name "intel" )
     elseif( ${my_cxx_compiler} MATCHES "pgCC")
        set( compiler_short_name "pgi" )
     endif()
     if( ${work_dir} MATCHES ".*[-]([0-9]+[.][0-9]+[.-][0-9]+).*" )
        string( REGEX REPLACE ".*[-]([0-9]+[.][0-9]+[.-][0-9]+).*" "\\1"
           compiler_version ${work_dir} )
     endif()
     if( NOT "${compiler_version}x" STREQUAL "x" )
        set( compiler_short_name "${compiler_short_name}-${compiler_version}" )
     endif()
  endif()

  # maybe just gcc?
  if( WIN32 )
    if( "$ENV{dirext}" MATCHES "x64" )
      set( CTEST_BUILD_NAME "Win64_${CTEST_BUILD_CONFIGURATION}" )
    else()
      set( CTEST_BUILD_NAME "Win32_${CTEST_BUILD_CONFIGURATION}" )
    endif()
  elseif( APPLE ) # OS/X
    set( CTEST_BUILD_NAME "OSX_${compiler_short_name}_${CTEST_BUILD_CONFIGURATION}" )
  else() # Unix
    set( CTEST_BUILD_NAME "Linux64_${compiler_short_name}_${CTEST_BUILD_CONFIGURATION}" )
  endif()

  # Default is no Coverage Analysis
  if( ${CTEST_SCRIPT_ARG} MATCHES Coverage )
    if( ${CTEST_BUILD_CONFIGURATION} MATCHES Release OR
        ${CTEST_BUILD_CONFIGURATION} MATCHES MinSizeRel )
      message( FATAL_ERROR "Cannot run coverage for \"Release\" mode builds." )
    endif()
    set( ENABLE_C_CODECOVERAGE ON )
    set( ENABLE_Fortran_CODECOVERAGE ON )
    set( CTEST_BUILD_NAME "${CTEST_BUILD_NAME}_Cov" )
  endif()
  
  # For Experimental builds, use launchers and parallel builds.
  if( ${CTEST_SCRIPT_ARG} MATCHES Experimental )
     if( UNIX )
        if( EXISTS "/proc/cpuinfo" )
           file( READ "/proc/cpuinfo" cpuinfo )
           # convert one big string into a set of strings, one per line
           string( REGEX REPLACE "\n" ";" cpuinfo ${cpuinfo} )
           set( proc_ids "" )
# consider using:
# include(ProcessorCount)
# ProcessorCount(DRACO_NUM_CORES)
           foreach( line ${cpuinfo} )
              if( ${line} MATCHES "processor" )
                 list( APPEND proc_ids ${line} )
              endif()
           endforeach()
           list( LENGTH proc_ids DRACO_NUM_CORES )
           set( MPIEXEC_MAX_NUMPROCS ${DRACO_NUM_CORES} CACHE STRING 
              "Number of cores on the local machine." )
        endif()
     endif()
  endif()

  if( ${drm_verbose} )    
    message("
CTEST_MODEL                 = ${CTEST_MODEL}
CTEST_BUILD_CONFIGURATION   = ${CTEST_BUILD_CONFIGURATION}
compiler_short_name         = ${compiler_short_name}
CTEST_BUILD_NAME            = ${CTEST_BUILD_NAME}
ENABLE_C_CODECOVERAGE       = ${ENABLE_C_CODECOVERAGE}
ENABLE_Fortran_CODECOVERAGE = ${ENABLE_Fortran_CODECOVERAGE}
CTEST_USE_LAUNCHERS         = ${CTEST_USE_LAUNCHERS}
MPIEXEC_MAX_NUMPROCS        = ${MPIEXEC_MAX_NUMPROCS}
")
  endif()
endmacro( parse_args )

# ------------------------------------------------------------
# Names of developer tools
# ------------------------------------------------------------
macro( find_tools )

  find_program( CTEST_CMD
    NAMES ctest
    HINTS
      "c:/Program Files (x86)/CMake 2.8/bin"
      # NO_DEFAULT_PATH
    )
  if( NOT EXISTS ${CTEST_CMD} )
    message( FATAL_ERROR "Cound not find ctest executable.(CTEST_CMD = ${CTEST_CMD})" )
  endif( NOT EXISTS ${CTEST_CMD} )

  find_program( CTEST_SVN_COMMAND
     NAMES svn
     HINTS
        "C:/Program Files (x86)/CollabNet Subversion"
        "C:/Program Files (x86)/CollabNet/Subversion Client"
        # NO_DEFAULT_PATH
     )
  set( CTEST_CVS_COMMAND ${CTEST_SVN_COMMAND} )
  if( NOT EXISTS "${CTEST_CVS_COMMAND}" )
    message( FATAL_ERROR "Cound not find cvs executable." )
  endif( NOT EXISTS "${CTEST_CVS_COMMAND}" )

  find_program( CTEST_CMAKE_COMMAND
    NAMES cmake
    HINTS
      "c:/Program Files (x86)/CMake 2.8/bin"
      # NO_DEFAULT_PATH
    )
  if( NOT EXISTS "${CTEST_CMAKE_COMMAND}" )
    message( FATAL_ERROR "Cound not find cmake executable." )
  endif( NOT EXISTS "${CTEST_CMAKE_COMMAND}" )

  find_program( MAKECOMMAND
    NAMES nmake make
    HINTS
      "C:/Program Files (x86)/Microsoft Visual Studio 11.0/VC/bin"
      "c:/Program Files (x86)/Microsoft Visual Studio 9.0/VC/bin"
      # NO_DEFAULT_PATH
    )
  if( NOT EXISTS "${MAKECOMMAND}" )
    message( FATAL_ERROR "Cound not find make/nmake executable." )
  endif()

  find_program( CTEST_MEMORYCHECK_COMMAND NAMES valgrind )
  set(          CTEST_MEMORYCHECK_COMMAND_OPTIONS  
     "-q --tool=memcheck --leak-check=full --trace-children=yes --suppressions=${CTEST_SCRIPT_DIRECTORY}/valgrind_suppress.txt --gen-suppressions=all" )
  # --show-reachable --num-callers=50
  # --suppressions=<filename>
  # --gen-suppressions=all|yes|no
  if( EXISTS ${CTEST_SCRIPT_DIRECTORY}/valgrind_suppress.txt )
     set( MEMORYCHECK_SUPPRESSIONS_FILE
        ${CTEST_SCRIPT_DIRECTORY}/valgrind_suppress.txt )
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
       execute_process(COMMAND ${COV01} -1
          RESULT_VARIABLE RES )
       if( RES )
          message(FATAL_ERROR "could not run cov01 -1")
       else()
          message("BullseyeCoverage turned on")
       endif()
    else()
       message( FATAL_ERROR 
          "Coverage requested, but bullseyecoverage's cov01 binary not in PATH."
          )
    endif()
  endif()

  # This breaks NMake Makefile builds because it is missing the -G"..."
  # set(CTEST_CONFIGURE_COMMAND "${CMAKE_COMMAND} \"${CTEST_SOURCE_DIRECTORY}\"")
  if( ${drm_verbose} )  
    message("
CTEST_CMD           = ${CTEST_CMD}
CTEST_CVS_COMMAND   = ${CTEST_CVS_COMMAND}
CTEST_CMAKE_COMMAND = ${CTEST_CMAKE_COMMAND}
MAKECOMMAND         = ${MAKECOMMAND}
CTEST_MEMORYCHECK_COMMAND         = ${CTEST_MEMORYCHECK_COMMAND}
MEMORYCHECK_SUPPRESSIONS_FILE     = ${MEMORYCHECK_SUPPRESSIONS_FILE}
CTEST_MEMORYCHECK_COMMAND_OPTIONS = ${CTEST_MEMORYCHECK_COMMAND_OPTIONS}
CTEST_CONFIGURE_COMMAND           = ${CTEST_CONFIGURE_COMMAND}

")
  endif()
endmacro( find_tools )

# ------------------------------------------------------------
# Setup regression steps:
# ------------------------------------------------------------
macro( setup_ctest_commands )
  # deprecated
endmacro( setup_ctest_commands )

# ------------------------------------------------------------
# Setup CVSROOT
# ------------------------------------------------------------
# macro( set_cvs_command projname )
#    if( EXISTS /ccs/codes/radtran/cvsroot )
#       set( CTEST_CVS_CHECKOUT
#          "${CTEST_CVS_COMMAND} -d /ccs/codes/radtran/cvsroot co -P -d source ${projname}" )
#    elseif( EXISTS /usr/projects/jayenne/regress/cvsroot )
#       set( CTEST_CVS_CHECKOUT
#          "${CTEST_CVS_COMMAND} -d /usr/projects/jayenne/regress/cvsroot co -P -d source ${projname}" )
#    else()
#       set( CTEST_CVS_CHECKOUT
#          "${CTEST_CVS_COMMAND} -d ccscs8:/ccs/codes/radtran/cvsroot co -P -d source ${projname}" )
#    endif()
# endmacro()
##---------------------------------------------------------------------------##
macro( set_svn_command svnpath )
   if( NOT EXISTS ${CTEST_SOURCE_DIRECTORY}/CMakeLists.txt )
      if( EXISTS /ccs/codes/radtran/svn )
         set( CTEST_CVS_CHECKOUT
            "${CTEST_CVS_COMMAND} checkout file:///ccs/codes/radtran/svn/${svnpath} source" )
         message("CTEST_CVS_CHECKOUT = ${CTEST_CVS_CHECKOUT}")
      elseif( EXISTS /usr/projects/jayenne/regress/svn )
         set( CTEST_CVS_CHECKOUT
            "${CTEST_CVS_COMMAND} checkout file:///usr/projects/jayenne/regress/svn/${svnpath} source" )
      else()
         set( CTEST_CVS_CHECKOUT
            "${CTEST_CVS_COMMAND} checkout svn+ssh://ccscs8/ccs/codes/radtran/svn/${svnpath} source" )
      endif()
   endif()
endmacro()

# ------------------------------------------------------------
# Setup for Code Coverage and LOC metrics
# ------------------------------------------------------------
macro( setup_for_code_coverage )
   if( "${sitename}" MATCHES "ccscs[89]" AND "$ENV{CXX}" MATCHES "g[+][+]" )
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
            execute_process(COMMAND "${COV01}" --on RESULT_VARIABLE RES)

            # Process and save lines of code 
            message( "Generating lines of code statistics...
 /home/regress/draco/regression/cloc
               --exclude-dir=heterogeneous,chimpy
               --exclude-list-file=/home/regress/draco/regression/cloc-exclude.cfg
               --exclude-lang=Text,Postscript
               --categorize=cloc-categorize.log 
               --counted=cloc-counted.log 
               --ignored=cloc-ignored.log 
               --progress-rate=0 
               --report-file=lines-of-code.log 
               --force-lang-def=/home/regress/draco/regression/cloc-lang.defs
               ${CTEST_SOURCE_DIRECTORY}
            ")
            execute_process( 
               COMMAND /home/regress/draco/regression/cloc
               --exclude-dir=heterogeneous,chimpy
               --exclude-list-file=/home/regress/draco/regression/cloc-exclude.cfg
               --exclude-lang=Text,Postscript
               --categorize=cloc-categorize.log 
               --counted=cloc-counted.log 
               --ignored=cloc-ignored.log 
               --progress-rate=0 
               --report-file=lines-of-code.log 
               --force-lang-def=/home/regress/draco/regression/cloc-lang.defs
               ${CTEST_SOURCE_DIRECTORY}
               #  --3 
               #  --diff
               WORKING_DIRECTORY ${CTEST_BINARY_DIRECTORY}
               )
            message( "Lines of code data at ${CTEST_BINARY_DIRECTORY}/lines-of-code.log")
            message( "Generating lines of code statistics (omitting test directories)
 /home/regress/draco/regression/cloc
               --exclude-dir=heterogeneous,chimpy,test
               --exclude-list-file=/home/regress/draco/regression/cloc-exclude.cfg
               --exclude-lang=Text,Postscript
               --categorize=cloc-categorize.log 
               --counted=cloc-counted.log 
               --ignored=cloc-ignored.log 
               --progress-rate=0 
               --report-file=lines-of-code-notest.log 
               --force-lang-def=/home/regress/draco/regression/cloc-lang.defs
               ${CTEST_SOURCE_DIRECTORY}
            ")
            execute_process( 
               COMMAND /home/regress/draco/regression/cloc
               --exclude-dir=heterogeneous,chimpy,test
               --exclude-list-file=/home/regress/draco/regression/cloc-exclude.cfg
               --exclude-lang=Text,Postscript
               --categorize=cloc-categorize.log 
               --counted=cloc-counted.log 
               --ignored=cloc-ignored.log 
               --progress-rate=0 
               --report-file=lines-of-code-notest.log 
               --force-lang-def=/home/regress/draco/regression/cloc-lang.defs
               ${CTEST_SOURCE_DIRECTORY}
               #  --3 
               #  --diff
               WORKING_DIRECTORY ${CTEST_BINARY_DIRECTORY}
               )
            message( "Lines of code data at ${CTEST_BINARY_DIRECTORY}/lines-of-code.log")

            set( CTEST_NOTES_FILES 
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
   if( "${sitename}" MATCHES "ccscs[89]" AND "$ENV{CXX}" MATCHES "g[+][+]" )
      if( ${CTEST_BUILD_CONFIGURATION} MATCHES Debug )
         if(ENABLE_C_CODECOVERAGE)
            message( "ctest_coverage( BUILD \"${CTEST_BINARY_DIRECTORY}\" )")
            ctest_coverage( BUILD "${CTEST_BINARY_DIRECTORY}" )
            execute_process(COMMAND "${COV01}" --off RESULT_VARIABLE RES)
         else()
            message( "ctest_memcheck( SCHEDULE_RANDOM ON )")
            ctest_memcheck(
               SCHEDULE_RANDOM ON 
               EXCLUDE_LABEL "nomemcheck")
         endif()
      endif()
   endif()
endmacro(process_cc_or_da)


# ------------------------------------------------------------
# Special default settings for a couple of platforms
# 
# ------------------------------------------------------------
macro(platform_customization)
   if( "${sitename}" MATCHES "Cielito" )
#      set( TOOLCHAIN_SETUP
#         "CMAKE_TOOLCHAIN_FILE:FILEPATH=/usr/projects/jayenne/regress/draco/config/Toolchain-catamount.cmake"
# )
      set(CT_CUSTOM_VARS 
"DRACO_LIBRARY_TYPE:STRING=STATIC
CMAKE_SYSTEM_NAME:STRING=Catamount
CMAKE_C_COMPILER:FILEPATH=cc
CMAKE_CXX_COMPILER:FILEPATH=CC 
CMAKE_Fortran_COMPILER:FILEPATH=ftn
MPIEXEC:FILEPATH=/usr/bin/aprun
MPIEXEC_NUMPROC_FLAG:STRING=-n
MPI_C_LIBRARIES:FILEPATH=
MPI_CXX_LIBRARIES:FILEPATH=
MPI_Fortran_LIBRARIES:FILEPATH=
MPI_C_INCLUDE_PATH:PATH=
MPI_CXX_INCLUDE_PATH:PATH=
MPI_Fortran_INCLUDE_PATH:PATH=")
   endif()
endmacro(platform_customization)

# ------------------------------------------------------------
# Special default settings for a couple of platforms
#
# Sets DRACO_DIR and CLUBIMC_DIR
# 
# ------------------------------------------------------------
macro(set_pkg_work_dir this_pkg dep_pkg)
   string( TOUPPER ${dep_pkg} dep_pkg_caps )
   # Assume that draco_work_dir is parallel to our current location.
   string( REPLACE ${this_pkg} ${dep_pkg} ${dep_pkg}_work_dir $ENV{work_dir} )
   if( "${dep_pkg}" MATCHES "draco" )
      string( REPLACE "cmake_jayenne/draco" "cmake_draco" 
         ${dep_pkg}_work_dir ${${dep_pkg}_work_dir} )      
   endif()
 
   find_file( ${dep_pkg}_target_dir
      NAMES README.${dep_pkg}
      HINTS
         # if DRACO_DIR is defined, use it.
         $ENV{DRACO_DIR}
         # regress account on ccscs8
         /home/regress/cmake_draco/${CTEST_MODEL}_${compiler_short_name}/${CTEST_BUILD_CONFIGURATION}/target
         # Try a path parallel to the work_dir
         ${${dep_pkg}_work_dir}/target
      )
   
   if( NOT EXISTS ${${dep_pkg}_target_dir} )
      message( FATAL_ERROR 
         "Could not locate the ${dep_pkg} installation directory. "
         "${dep_pkg}_target_dir = ${${dep_pkg}_target_dir}" )
   endif()
   get_filename_component( ${dep_pkg_caps}_DIR ${${dep_pkg}_target_dir} PATH )
   
   if( ENABLE_C_CODECOVERAGE )
      string( REPLACE ${CTEST_BUILD_CONFIGURATION} "Coverage" 
         ${dep_pkg_caps}_DIR ${${dep_pkg_caps}_DIR} )
   endif()
endmacro(set_pkg_work_dir)

