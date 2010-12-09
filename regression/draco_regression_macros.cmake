# Call this script from regress/Draco_*.cmake

# Sample script:
#----------------------------------------
# include( "${CTEST_SCRIPT_DIRECTORY}/draco_regression_macros.cmake" )
# set_defaults()
# parse_args()
# find_tools()
# setup_ctest_commands()
# set( CTEST_INITIAL_CACHE "
#VERBOSE:BOOL=ON
#BUILD_TESTS:BOOL=ON

#BUILDNAME:STRING=${build_name}
#CMAKE_BUILD_TYPE:STRING=${build_type}
#CMAKE_GENERATOR:STRING=${CMAKE_GENERATOR}
#CMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
#CMAKE_MAKE_PROGRAM:FILEPATH=${MAKECOMMAND}
#CVSCOMMAND:FILEPATH=${CTEST_CVS_COMMAND}
#SVNCOMMAND:FILEPATH=${CTEST_CVS_COMMAND}
#MAKECOMMAND:FILEPATH=${MAKECOMMAND}
#SITE:STRING=${sitename}
#VENDOR_DIR:PATH=${VENDOR_DIR}
#${CTEST_INITIAL_CACHE_EXTRAS}
#")
#set( CTEST_ENVIRONMENT
#  FC=/opt/pathscale/bin/pathf90
#  CXX=
#  CC=
#  VERBOSE=ON
#  CTEST_OUTPUT_ON_FAILURE=ON
#)
#----------------------------------------

# available commands for ctest:
# exec_program
# execute_proces
# file
# find_file
# find_library
# find_package
# find_path
# find_program
# macro
# site_name
# unset
# string

# ------------------------------------------------------------
# Defaults (override with optional arguments)
# ------------------------------------------------------------
macro( set_defaults )

  # parse arguments
  if( ${ARGV} MATCHES QUIET )
    set( quiet_mode ON )
  endif( ${ARGV} MATCHES QUIET )

  # Prerequisits
  if( NOT work_dir )
    if( EXISTS "$ENV{work_dir}" )
      set( work_dir $ENV{work_dir} )
    else( EXISTS "$ENV{work_dir}" )
      if( EXISTS "$ENV{PWD}" )
        set( work_dir $ENV{PWD} )
        message( "
Warning:  work_dir has been set to pwd.
Set work_dir in your environment if you want to use a different
location.
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
  set( CTEST_SOURCE_DIRECTORY "${work_dir}/source" )
  set( CTEST_BINARY_DIRECTORY "${work_dir}/build"  )
  set( CMAKE_INSTALL_PREFIX   "${work_dir}/target"  )
  
  # Default is "Experimental." Special builds are "Nightly" or "Continuous"
  set( dashboard_type "Experimental" ) 

  # Default is "Release." 
  # Special types are "Debug," "RelWithDebInfo" or "MinSizeRel"
  set( build_type "Release" )
  
  # should ctest wipe the binary tree before running
  set( CTEST_START_WITH_EMPTY_BINARY_DIRECTORY TRUE )

  if( EXISTS "$ENV{VENDOR_DIR}" )
    set(VENDOR_DIR $ENV{VENDOR_DIR})
  endif()
  find_path( VENDOR_DIR
    ChangeLog
    PATHS
      /ccs/codes/radtran/vendors/Linux64
      /radiative/vendors/${CMAKE_SYSTEM_PROCESSOR}-${CMAKE_SYSTEM_NAME}
      c:/vendors/${CMAKE_SYSTEM_PROCESSOR}-${CMAKE_SYSTEM_NAME}
      c:/vendors
      )
#      /radiative/vendors/${CMAKE_SYSTEM_PROCESSOR}-${CMAKE_SYSTEM_NAME}
  
  if( NOT EXISTS ${VENDOR_DIR} )
    message( FATAL_ERROR "VENDOR_DIR = ${VENDOR_DIR} was not found." )
  endif( NOT EXISTS ${VENDOR_DIR} )

  if( WIN32 )
    # add option for "NMake Makefiles JOM"?
    set( CMAKE_GENERATOR "NMake Makefiles" )
  else( WIN32 )
    set( CMAKE_GENERATOR "Unix Makefiles" )
  endif( WIN32 )      

  # Set the sitename, but strip any domain information
  site_name( sitename )
  string( REGEX REPLACE "([A-z0-9]+).*" "\\1" sitename ${sitename} )
  # Treat all Infinitron nodes as infinitron
  string( REGEX REPLACE "n00[0-9]" "infinitron" sitename ${sitename} )

  set( ENABLE_C_CODECOVERAGE OFF )
  set( ENABLE_Fortran_CODECOVERAGE OFF )

  if( NOT quiet_mode )
    message("
sitename = ${sitename}

CTEST_SCRIPT_NAME      = ${CTEST_SCRIPT_NAME}
CTEST_SCRIPT_DIRECTORY = ${CTEST_SCRIPT_DIRECTORY}
CTEST_SOURCE_DIRECTORY = ${CTEST_SOURCE_DIRECTORY}
CTEST_BINARY_DIRECTORY = ${CTEST_BINARY_DIRECTORY}
CMAKE_INSTALL_PREFIX   = ${CMAKE_INSTALL_PREFIX}
CMAKE_GENERATOR        = ${CMAKE_GENERATOR}
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
    set( dashboard_type "Nightly" )
  elseif( ${CTEST_SCRIPT_ARG} MATCHES Continuous )
    set( dashboard_type "Continuous" )
  endif( ${CTEST_SCRIPT_ARG} MATCHES Nightly )
  
  # Default is "Release." 
  # Special types are "Debug," "RelWithDebInfo" or "MinSizeRel"
  if( ${CTEST_SCRIPT_ARG} MATCHES Debug )
    set( build_type "Debug" )
  elseif( ${CTEST_SCRIPT_ARG} MATCHES RelWithDebInfo )
    set( build_type "RelWithDebInfo" )
  elseif( ${CTEST_SCRIPT_ARG} MATCHES MinSizeRel )
    set( build_type "MinSizeRel" )
  endif( ${CTEST_SCRIPT_ARG} MATCHES Debug )
  
  # maybe just gcc?
  if( WIN32 )
    if( "$ENV{dirext}" MATCHES "x64" )
      set( build_name "Win64_${build_type}" )
    else()
      set( build_name "Win32_${build_type}" )
    endif()
  else() # Unix
    set( build_name "Linux64_${build_type}" )
  endif()

  # Default is no Coverage Analysis
  if( ${CTEST_SCRIPT_ARG} MATCHES Coverage )
    if( ${build_type} MATCHES Release OR
        ${build_type} MATCHES MinSizeRel )
      message( FATAL_ERROR "Cannot run coverage for \"Release\" mode builds." )
    endif()
    set( enable_coverage ON )
    set( build_name "${build_name}_Cov" )
  endif()
  
  if( NOT quiet_mode )
    message("
dashboard_type = ${dashboard_type}
build_type     = ${build_type}
build_name     = ${build_name}
enable_coverage= ${enable_coverage}
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
    message( FATAL_ERROR "Cound not find ctest executable." )
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
  endif( NOT EXISTS "${MAKECOMMAND}" )

  find_program( CTEST_MEMORYCHECK_COMMAND NAMES valgrind )
  find_program( MEMORYCHECK_COMMAND NAMES valgrind )

  set( MEMORYCHECK_COMMAND_OPTIONS "--leak-check=yes" CACHE STRING 
  "Options for memorycheck tool (valgrind)." )
# "--leak-check=yes --num-callers=8 --show-reachable=yes"
# "--leak-check=full"


  if( NOT quiet_mode )
    message("
CTEST_CMD           = ${CTEST_CMD}
CTEST_CVS_COMMAND   = ${CTEST_CVS_COMMAND}
CTEST_CMAKE_COMMAND = ${CTEST_CMAKE_COMMAND}
MAKECOMMAND         = ${MAKECOMMAND}
")
  endif( NOT quiet_mode )

endmacro( find_tools )

# ------------------------------------------------------------
# Setup regression steps:
# ------------------------------------------------------------
macro( setup_ctest_commands )
  # which ctest command to use for running the dashboard
  # Do the normal Start,Upate,Configure,Build,Test,Submit

  # ${dashboard_type}Start
  # ${dashboard_type}Update
  # ${dashboard_type}Configure
  # ${dashboard_type}Build
  # ${dashboard_type}Test
  # ${dashboard_type}Submit
  # Skip: 
  # ${dashboard_type}Coverage
  # ${dashboard_type}MemCheck

    # parse arguments
  if( ${ARGV} MATCHES QUIET )
    set( quiet_mode ON )
  endif( ${ARGV} MATCHES QUIET )

  set( CTEST_COMMAND
    "${CTEST_CMD} -D ${dashboard_type}Start -D ${dashboard_type}Update -D ${dashboard_type}Configure -D ${dashboard_type}Build -D ${dashboard_type}Test -D ${dashboard_type}Submit"
    )
  
  if( ${dashboard_type} MATCHES Continuous )
    # only run tests 1 and 2.
    set( CTEST_COMMAND "${CTEST_CMD} -D ${dashboard_type} -I1,2" )
    # Only recompile minimum set of files after first build.
    set( CTEST_START_WITH_EMPTY_BINARY_DIRECTORY_ONCE TRUE )
    set( CTEST_START_WITH_EMPTY_BINARY_DIRECTORY FALSE )
    # How long are the continuous builds active? (only used for
    # Continuous mode.)
    set( CTEST_CONTINUOUS_DURATION 180 ) # minutes
    # How long to wait before starting next build? (only used for
    # Continuous mode.)
    set( CTEST_CONTINUOUS_MINIMUM_INTERVAL 30 ) # minutes
  endif( ${dashboard_type} MATCHES Continuous )

  # Do Dynamic Analysis (MemCheck) when Debug, but not Coverage.
  if( UNIX )
    if( ${build_type} MATCHES Debug )
      if( enable_coverage )
        set( CTEST_COMMAND
          "${CTEST_COMMAND} -D ${dashboard_type}Coverage -D ${dashboard_type}Submit" )
        set( ENABLE_C_CODECOVERAGE ON )
        set( ENABLE_Fortran_CODECOVERAGE ON )
#        set( CTEST_INITIAL_CACHE_EXTRAS "ENABLE_C_CODECOVERAGE:BOOL=ON" )

      else( enable_coverage )
        set( CTEST_COMMAND
          "${CTEST_COMMAND} -D ${dashboard_type}MemCheck -D ${dashboard_type}Submit"
          )
      endif( enable_coverage )
    endif( ${build_type} MATCHES Debug )
  endif( UNIX )

  if( NOT quiet_mode )
    message("

CTEST_COMMAND = ${CTEST_COMMAND}
")
  endif( NOT quiet_mode )

endmacro( setup_ctest_commands )

