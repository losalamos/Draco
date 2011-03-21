#-----------------------------*-cmake-*----------------------------------------#
# file   draco_regression_macros.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 Nov 22
# brief  Helper macros for setting up a CTest/CDash regression system
# note   Copyright © 2010-2011 Los Alamos National Security
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

# Ref: http://www.cmake.org/Wiki/CMake_Testing_With_CTest
#      http://www.cmake.org/Wiki/CMake_Scripting_of_CTest


# Call this script from regress/Draco_*.cmake

# Sample script:
#----------------------------------------
# set( CTEST_PROJECT_NAME "Draco" )
# include( "${CTEST_SCRIPT_DIRECTORY}/draco_regression_macros.cmake" )
# set_defaults()
# parse_args()
# find_tools()
# set( CTEST_INITIAL_CACHE "
# CMAKE_VERBOSE_MAKEFILE:BOOL=ON
# CMAKE_BUILD_TYPE:STRING=${CTEST_BUILD_CONFIGURATION}
# CMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
# CTEST_CMAKE_GENERATOR:STRING=${CTEST_CMAKE_GENERATOR}
# CTEST_USE_LAUNCHERS:STRING=${CTEST_USE_LAUNCHERS}
# ENABLE_C_CODECOVERAGE:BOOL=${ENABLE_C_CODECOVERAGE}
# ENABLE_Fortran_CODECOVERAGE:BOOL=${ENABLE_Fortran_CODECOVERAGE}
# VENDOR_DIR:PATH=/ccs/codes/radtran/vendors/Linux64
# ")
# ctest_empty_binary_directory( ${CTEST_BINARY_DIRECTORY} )
# file( WRITE ${CTEST_BINARY_DIRECTORY}/CMakeCache.txt ${CTEST_INITIAL_CACHE} )
# set( VERBOSE ON )
# set( CTEST_OUTPUT_ON_FAILURE ON )
# setup_ctest_commands()
#execute_process( 
#   COMMAND           ${CMAKE_MAKE_PROGRAM} install
#   WORKING_DIRECTORY ${CTEST_BINARY_DIRECTORY}
#   )
#----------------------------------------

# ------------------------------------------------------------
# Defaults (override with optional arguments)
# ------------------------------------------------------------
macro( set_defaults )

  # parse arguments
  if( ${ARGV} MATCHES QUIET )
    set( quiet_mode ON )
  endif( ${ARGV} MATCHES QUIET )

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
  if( ${sitename} MATCHES "yr[a-d]+[0-9]+" )
     set( sitename "YellowRail" )
  elseif( ${sitename} MATCHES "tu[a-d]+[0-9]+" )
     set( sitename "Turing" )
  endif()
  # string( REGEX REPLACE "n00[0-9]" "infinitron" sitename ${sitename} )
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
  else( WIN32 )
    set( CTEST_CMAKE_GENERATOR "Unix Makefiles" )
  endif( WIN32 )      

  #Only works for makefile generators, but gives us pretty output.
#  if( ${CTEST_CMAKE_GENERATOR} STREQUAL "Unix Makefiles" )
#     set( CTEST_USE_LAUNCHERS 1 )
#  else()
     set( CTEST_USE_LAUNCHERS 0 )
#  endif()

  set( ENABLE_C_CODECOVERAGE OFF )
  set( ENABLE_Fortran_CODECOVERAGE OFF )

  # Dashboard setup (in place of CTestConfig.cmake)
  if( NOT CTEST_PROJECT_NAME )
     set( CTEST_PROJECT_NAME "UnknownProject")
  endif()
  set( CTEST_NIGHTLY_START_TIME "00:10:00 MST")
  
  set( CTEST_DROP_METHOD "http")
  set( CTEST_DROP_SITE "coder.lanl.gov")
  set( CTEST_DROP_LOCATION 
     "/cdash/submit.php?project=${CTEST_PROJECT_NAME}" )
  set( CTEST_DROP_SITE_CDASH TRUE )
  set( CTEST_CURL_OPTIONS CURLOPT_SSL_VERIFYPEER_OFF )

# MCATK settings
#  set( CTEST_DROP_METHOD "https")
#  set( CTEST_DROP_SITE "cdash.lanl.gov")
#  set( CTEST_DROP_LOCATION "/submit.php?project=test" )
#  set( CTEST_DROP_SITE_CDASH TRUE )
#  set( CTEST_CURL_OPTIONS CURLOPT_SSL_VERIFYPEER_OFF )

  # if( UNIX )
  #    if( EXISTS "/proc/cpuinfo" )
  #       file( READ "/proc/cpuinfo" cpuinfo )
  #       # convert one big string into a set of strings, one per line
  #       string( REGEX REPLACE "\n" ";" cpuinfo ${cpuinfo} )
  #       set( proc_ids "" )
  #       foreach( line ${cpuinfo} )
  #          if( ${line} MATCHES "processor" )
  #             list( APPEND proc_ids ${line} )
  #          endif()
  #       endforeach()
  #       list( LENGTH proc_ids DRACO_NUM_CORES )
  #       set( MPIEXEC_MAX_NUMPROCS ${DRACO_NUM_CORES} CACHE STRING 
  #          "Number of cores on the local machine." )
  #    endif()
  # endif()
  set( MPIEXEC_MAX_NUMPROCS 1 CACHE STRING  "Number of cores on the local machine." )

  if( EXISTS "$ENV{VENDOR_DIR}" )
    set(VENDOR_DIR $ENV{VENDOR_DIR})
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

   set( VERBOSE ON )
   set( CTEST_OUTPUT_ON_FAILURE ON )

# Consider setting the following:

# SET(CTEST_CUSTOM_ERROR_PRE_CONTEXT 20)
# SET(CTEST_CUSTOM_ERROR_POST_CONTEXT 20)
# SET(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_ERRORS 100)
# SET(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_WARNINGS 100)
#     CTEST_START_WITH_EMPTY_BINARY_DIRECTORY
#     CTEST_CONTINUOUS_DURATION
#     CTEST_CONTINUOUS_MINIMUM_INTERVAL
#find_program(CTEST_GIT_COMMAND NAMES git PATHS "C:\\Program Files\\Git\\bin")
#find_program(CTEST_COVERAGE_COMMAND NAMES gcov)

  # Echo settings
  
  if( NOT quiet_mode )
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
CTEST_USE_LAUNCHERS    = ${CTEST_USE_LAUNCHERS}
ENABLE_C_CODECOVERAGE  = ${ENABLE_C_CODECOVERAGE}
ENABLE_Fortran_CODECOVERAGE = ${ENABLE_Fortran_CODECOVERAGE}
CTEST_NIGHTLY_START_TIME  = ${CTEST_NIGHTLY_START_TIME}
CTEST_DROP_METHOD         = ${CTEST_DROP_METHOD}
CTEST_DROP_SITE           = ${CTEST_DROP_SITE}
CTEST_DROP_LOCATION       = ${CTEST_DROP_LOCATION}
CTEST_DROP_SITE_CDASH     = ${CTEST_DROP_SITE_CDASH}
CTEST_CURL_OPTIONS        = ${CTEST_CURL_OPTIONS}
MPIEXEC_MAX_NUMPROCS      = ${MPIEXEC_MAX_NUMPROCS}
VENDOR_DIR                = ${VENDOR_DIR}
")
  endif( NOT quiet_mode )

endmacro( set_defaults )

# ------------------------------------------------------------
# Parse Arguments
# ------------------------------------------------------------
macro( parse_args )

  # parse arguments
  if( ${ARGV} MATCHES QUIET )
    set( quiet_mode ON )
  endif( ${ARGV} MATCHES QUIET )

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
  
  set( compiler_short_name "gcc" )
  if( $ENV{CXX} MATCHES "pgCC" )
     set( compiler_short_name "pgi" )
  elseif($ENV{CXX} MATCHES "icpc" )
     set( compiler_short_name "intel" )
  endif()

  # maybe just gcc?
  if( WIN32 )
    if( "$ENV{dirext}" MATCHES "x64" )
      set( CTEST_BUILD_NAME "Win64_${CTEST_BUILD_CONFIGURATION}" )
    else()
      set( CTEST_BUILD_NAME "Win32_${CTEST_BUILD_CONFIGURATION}" )
    endif()
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
  
  if( NOT quiet_mode )
    message("
CTEST_MODEL               = ${CTEST_MODEL}
CTEST_BUILD_CONFIGURATION = ${CTEST_BUILD_CONFIGURATION}
compiler_short_name       = ${compiler_short_name}
CTEST_BUILD_NAME          = ${CTEST_BUILD_NAME}
ENABLE_C_CODECOVERAGE     = ${ENABLE_C_CODECOVERAGE}
ENABLE_Fortran_CODECOVERAGE = ${ENABLE_Fortran_CODECOVERAGE}
")
  endif()

endmacro( parse_args )

# ------------------------------------------------------------
# Names of developer tools
# ------------------------------------------------------------
macro( find_tools )

  # parse arguments
  if( ${ARGV} MATCHES QUIET )
    set( quiet_mode ON )
  endif( ${ARGV} MATCHES QUIET )

  find_program( CTEST_CMD
    NAMES ctest
    HINTS
      "c:/Program Files (x86)/CMake 2.8/bin"
      # NO_DEFAULT_PATH
    )
  if( NOT EXISTS ${CTEST_CMD} )
    message( FATAL_ERROR "Cound not find ctest executable.(CTEST_CMD = ${CTEST_CMD})" )
  endif( NOT EXISTS ${CTEST_CMD} )

  find_program( CTEST_CVS_COMMAND
    NAMES cvs
    HINTS
      "C:/Program Files (x86)/CollabNet Subversion"
      "C:/Program Files (x86)/CollabNet/Subversion Client"
      # NO_DEFAULT_PATH
    )
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
      "c:/Program Files (x86)/Microsoft Visual Studio 9.0/VC/bin"
      # NO_DEFAULT_PATH
    )
  if( NOT EXISTS "${MAKECOMMAND}" )
    message( FATAL_ERROR "Cound not find make/nmake executable." )
  endif()

  find_program( CTEST_MEMORYCHECK_COMMAND NAMES valgrind )
  set(          CTEST_MEMORYCHECK_COMMAND_OPTIONS  
     "-q --tool=memcheck --leak-check=full --trace-children=yes --error-limit=100 --suppressions=${CTEST_SCRIPT_DIRECTORY}/valgrind_suppress.txt" )
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
          message(STATUS "BullseyeCoverage turned on")
       endif()
    else()
       message( FATAL_ERROR 
          "Coverage requested, but bullseyecoverage's cov01 binary not in PATH."
          )
    endif()
  endif()

#   if( MPIEXEC_MAX_NUMPROCS )
#      set( CMAKE_BUILD_COMMAND "gmake -j ${MPIEXEC_MAX_NUMPROCS}" ) # install | check
# #     set( CTEST_CMD "${CTEST_CMD} -j ${MPIEXEC_MAX_NUMPROCS}" )
#      set( MAKECOMMAND "${MAKECOMMAND} -j ${MPIEXEC_MAX_NUMPROCS}" )
# #    set( CTEST_BUILD_COMMAND ${CMAKE_BUILD_COMMAND} )     
#   endif()

  set(CTEST_CONFIGURE_COMMAND "${CMAKE_COMMAND} \"${CTEST_SOURCE_DIRECTORY}\"")

  if( NOT quiet_mode )
    message("
CTEST_CMD           = ${CTEST_CMD}
CTEST_CVS_COMMAND   = ${CTEST_CVS_COMMAND}
CTEST_CMAKE_COMMAND = ${CTEST_CMAKE_COMMAND}
MAKECOMMAND         = ${MAKECOMMAND}
CTEST_MEMORYCHECK_COMMAND     = ${CTEST_MEMORYCHECK_COMMAND}
MEMORYCHECK_SUPPRESSIONS_FILE = ${MEMORYCHECK_SUPPRESSIONS_FILE}
CTEST_MEMORYCHECK_COMMAND_OPTIONS = ${CTEST_MEMORYCHECK_COMMAND_OPTIONS}
beyedir                       = ${beyedir}
CTEST_CONFIGURE_COMMAND       = ${CTEST_CONFIGURE_COMMAND}
")
  endif( NOT quiet_mode )

endmacro( find_tools )

# ------------------------------------------------------------
# Setup regression steps:
# ------------------------------------------------------------
macro( setup_ctest_commands )
  # which ctest command to use for running the dashboard
  # Do the normal Start,Upate,Configure,Build,Test,Submit

    # parse arguments
  if( ${ARGV} MATCHES QUIET )
    set( quiet_mode ON )
  endif( ${ARGV} MATCHES QUIET )

  #
  # Drive the problem (www.cmake.org/cmake/help/ctest-2-8=docs.html)
  #

  message(STATUS "ctest_start( ${CTEST_MODEL} )")
#  ctest_start( ${CTEST_MODEL} )
#  message(STATUS  "ctest_update()"  )
#  ctest_update()
#   if( "$ENV{CXX}" MATCHES "g[+][+]" )
#      if( ${CTEST_BUILD_CONFIGURATION} MATCHES Debug )
#         if(ENABLE_C_CODECOVERAGE)
#            configure_file( 
#               ${CTEST_SCRIPT_DIRECTORY}/covclass_cmake.cfg
#               ${CTEST_BINARY_DIRECTORY}/covclass_cmake.cfg 
#               @ONLY )
#            set( ENV{COVDIRCFG}   ${CTEST_BINARY_DIRECTORY}/covclass_cmake.cfg )
#            set( ENV{COVFNCFG}    ${CTEST_BINARY_DIRECTORY}/covclass_cmake.cfg )
#            set( ENV{COVCLASSCFG} ${CTEST_BINARY_DIRECTORY}/covclass_cmake.cfg )
#            set( ENV{COVSRCCFG}   ${CTEST_BINARY_DIRECTORY}/covclass_cmake.cfg )
#            set( ENV{COVFILE}     ${CTEST_BINARY_DIRECTORY}/CMake.cov )
#            execute_process(COMMAND "${COV01}" --on
#               RESULT_VARIABLE RES)
#         endif()
#      endif()
#   endif()
#   message(STATUS "ctest_configure()" )
#   ctest_configure() # LABELS label1 [label2]
#   message(STATUS "ctest_build()" )
#   ctest_build()
#   message(STATUS "ctest_test(PARALLEL_LEVEL ${MPIEXEC_MAX_NUMPROCS} SCHEDULE_RANDOM ON )" )
#   ctest_test( 
#      PARALLEL_LEVEL ${MPIEXEC_MAX_NUMPROCS} 
#      SCHEDULE_RANDOM ON ) 
# #     EXCLUDE_LABEL "nomemcheck" )

#   if( "$ENV{CXX}" MATCHES "g[+][+]" )
#      if( ${CTEST_BUILD_CONFIGURATION} MATCHES Debug )
#         if(ENABLE_C_CODECOVERAGE)
#            message(STATUS "ctest_coverage( BUILD \"${CTEST_BINARY_DIRECTORY}\" )")
#            ctest_coverage( BUILD "${CTEST_BINARY_DIRECTORY}" )  # LABLES "scalar tests" 
#            execute_process(COMMAND "${COV01}" --off
#               RESULT_VARIABLE RES)
#         else()
#            message(STATUS "ctest_memcheck( PARALLEL_LEVEL ${MPIEXEC_MAX_NUMPROCS} SCHEDULE_RANDOM ON )")
#            ctest_memcheck(
#               SCHEDULE_RANDOM ON 
#               EXCLUDE_LABEL "nomemcheck")
# #              PARALLEL_LEVEL  ${MPIEXEC_MAX_NUMPROCS} 
#         endif()
#      endif()
#   endif()
#   message(STATUS "ctest_submit()")
#   ctest_submit()
     
endmacro( setup_ctest_commands )
