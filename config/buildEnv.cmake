#-----------------------------*-cmake-*----------------------------------------#
# file   buildEnv.cmake
# author Kelly Thompson <kgt@lanl.gov>
# date   2010 June 5
# brief  Default CMake build parameters
# note   Copyright (C) 2010-2012 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Build Parameters
#------------------------------------------------------------------------------#
macro( dbsSetDefaults )

  # if undefined, force build_type to "release"
  if( "${CMAKE_BUILD_TYPE}x" STREQUAL "x" )
    set( CMAKE_BUILD_TYPE "Debug" CACHE STRING 
       "Release, Debug, RelWithDebInfo" FORCE )
  endif( "${CMAKE_BUILD_TYPE}x" STREQUAL "x" )

  # constrain pull down values in cmake-gui
  set_property( CACHE CMAKE_BUILD_TYPE 
     PROPERTY 
       STRINGS Release Debug MinSizeRel RelWithDebInfo )

  # Provide default value for install_prefix
  if( "${CMAKE_INSTALL_PREFIX}" STREQUAL "/usr/local" OR
      "${CMAKE_INSTALL_PREFIX}" MATCHES "C:/Program Files" )
    set( CMAKE_INSTALL_PREFIX "${Draco_BINARY_DIR}/install" CACHE PATH 
       "Install path prefix, prepended onto install directories" FORCE)
  endif()
  
  mark_as_advanced( EXECUTABLE_OUTPUT_PATH )
  mark_as_advanced( LIBRARY_OUTPUT_PATH )
  mark_as_advanced( DART_TESTING_TIMEOUT )
  
  # Establish the sitename
  site_name( current_site )
  set( SITE ${current_site} CACHE INTERNAL 
     "Name of the computer/site where compile is being run." FORCE )
  mark_as_advanced( current_site )

  # Option for solution folders for GUI-based development environments
  if( ${CMAKE_GENERATOR} MATCHES "Visual Studio" )
     option( BUILD_USE_SOLUTION_FOLDERS "Enable grouping of projects in VS" ON )
     set_property( GLOBAL PROPERTY USE_FOLDERS ${BUILD_USE_SOLUTION_FOLDERS} )
  endif()

  # Special setup for Visual Studio builds.
  if(WIN32 AND NOT UNIX AND NOT BORLAND AND NOT MINGW )
     set( CMAKE_SUPPRESS_REGENERATION ON )  
  endif()
  
  # Design-by-Contract

  #   Insist() assertions    : always on
  #   Require() preconditions: add +1 to DBC_LEVEL
  #   Check() assertions     : add +2 to DBC_LEVEL
  #   Ensure() postconditions: add +4 to DBC_LEVEL
  #   Do not throw on error  : add +8 to DBC_LEVEL
  string( TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE )
  if( "${CMAKE_BUILD_TYPE}" MATCHES "RELEASE" )
     set( DRACO_DBC_LEVEL "0" CACHE STRING "Design-by-Contract (0-15)?" )
  else()
     set( DRACO_DBC_LEVEL "7" CACHE STRING "Design-by-Contract (0-15)?" )
  endif()
  # provide a constrained drop down menu in cmake-gui
  set_property( CACHE DRACO_DBC_LEVEL PROPERTY STRINGS 
     0 1 2 3 4 5 6 7 9 10 11 12 13 14 15)
  
  if( CMAKE_CONFIGURATION_TYPES )
     set( NESTED_CONFIG_TYPE -C "${CMAKE_CFG_INTDIR}" )
  else()
     if( CMAKE_BUILD_TYPE )
        set( NESTED_CONFIG_TYPE -C "${CMAKE_BUILD_TYPE}" )
     else()
        set( NESTED_CONFIG_TYPE )
     endif()
  endif()
  
  # Enable parallel build for Eclipse:
  set( CMAKE_ECLIPSE_MAKE_ARGUMENTS "-j ${MPIEXEC_MAX_NUMPROCS}" )

endmacro()

##---------------------------------------------------------------------------##
## dbsInitExportTargets
##
## These fields are constructed during Draco configure and are
## saved/installed to lib/cmake/draco-X.X/draco-config.cmake.
##---------------------------------------------------------------------------##
macro( dbsInitExportTargets PREFIX )
  # Data for exporting during install
   set( ${PREFIX}_LIBRARIES ""  CACHE INTERNAL "List of draco targets" FORCE)
   set( ${PREFIX}_PACKAGE_LIST ""  CACHE INTERNAL
      "List of known package targets" FORCE)
   set( ${PREFIX}_TPL_LIST ""  CACHE INTERNAL 
      "List of third party libraries known by this package" FORCE)
   set( ${PREFIX}_TPL_INCLUDE_DIRS ""  CACHE
      INTERNAL "List of include paths used by this package to find thrid party vendor header files." 
      FORCE)
   set( ${PREFIX}_TPL_LIBRARIES ""  CACHE INTERNAL
      "List of third party libraries used by this package." FORCE)
endmacro()

#------------------------------------------------------------------------------#
# Save some build parameters for later use by --config options
#------------------------------------------------------------------------------#
macro( dbsConfigInfo )

   # create a generator to capture DBS_BUILD_DATE
   # set( DBS_CONFIGURE_DATE "unknown" )

   # Do we have access to the all the required system headers?
   #include( CheckIncludeFiles )
   #check_include_files( stdbool.h DBS_HAVE_STDBOOL_H )
   #check_include_files( stdint.h TSP_HAVE_STDINT_H )
   #check_include_files( regex.h TSP_HAVE_REGEX_H )

   # Do we have some special symbols
   #include( CheckSymbolExists )
   #check_symbol_exists( cbrt math.h DBS_HAVE_CBRT )
   
   # We need this string in all caps to capture the compiler flags correctly.
   string( TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE )
   
   set( DBS_F90_COMPILER "${CMAKE_Fortran_COMPILER}" )
   set( DBS_F90_COMPILER_VER "${CMAKE_Fortran_COMPILER_VERSION}" )
   set( DBS_F90_OPT_LEVEL "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}}" )

   set( DBS_C_COMPILER "${CMAKE_C_COMPILER}" )
   set( DBS_C_OPT_LEVEL "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}}" )
   if( NOT DBS_C_COMPILER_VER )
      set( DBS_C_COMPILER_VER " " )
   endif( NOT DBS_C_COMPILER_VER )
   
   set( DBS_CXX_COMPILER "${CMAKE_CXX_COMPILER}" )
   set( DBS_CXX_OPT_LEVEL "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}")
   if( NOT DBS_CXX_COMPILER_VER )
      set( DBS_CXX_COMPILER_VER " ")
   endif( NOT DBS_CXX_COMPILER_VER )

   set( DBS_OPERATING_SYSTEM "${CMAKE_SYSTEM_NAME}")
   set( DBS_OPERATING_SYSTEM_VER "${CMAKE_SYSTEM}")

   # Suppliment with system commands as needed:
   if( UNIX )

      # Get some extra version information if this is RedHat.
      if( EXISTS "/etc/redhat-release" )
         file( READ "/etc/redhat-release" redhat_version )
         string( STRIP "${redhat_version}" redhat_version )
         set( DBS_OPERATING_SYSTEM_VER "${redhat_version} (${CMAKE_SYSTEM})" )
      endif( EXISTS "/etc/redhat-release" )

      # How many local cores
      if( EXISTS "/proc/cpuinfo" )
         file( READ "/proc/cpuinfo" cpuinfo )
         # string( STRIP "${cpuinfo}" cpuinfo )
         # convert one big string into a set of strings, one per line
         string( REGEX REPLACE "\n" ";" cpuinfo ${cpuinfo} )
         set( proc_ids "" )
         foreach( line ${cpuinfo} )
            if( ${line} MATCHES "processor" )
               list( APPEND proc_ids ${line} )
            endif()
         endforeach()
         list( LENGTH proc_ids DRACO_NUM_CORES )
         set( MPIEXEC_MAX_NUMPROCS ${DRACO_NUM_CORES} CACHE STRING 
            "Number of cores on the local machine." )
      endif()
      
   elseif() # WIN32

      # OS version information
      # Windows XP
      GET_FILENAME_COMPONENT( win_prod_name "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\Windows NT\\CurrentVersion;ProductName]" NAME )
      GET_FILENAME_COMPONENT( win_sp "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\Windows NT\\CurrentVersion;CSDVersion]" NAME )
      get_filename_component( win_ver "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\Windows NT\\CurrentVersion;CurrentVersion]" NAME )
      get_filename_component( win_build "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\Windows NT\\CurrentVersion;CurrentBuildNumber]" NAME )
      get_filename_component( win_buildlab "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\Windows NT\\CurrentVersion;BuildLab]" NAME )

      set( DBS_OPERATING_SYSTEM "${win_prod_name}" )
      set( DBS_OPERATING_SYSTEM_VER "${win_ver}.${win_build} (${win_buildlab})" )
      if( EXISTS "C:/Program Files (x86)" )
         set( windows_bits "x64" )
         set( DBS_OPERATING_SYSTEM "${win_prod_name} (${windows_bits}) ${win_sp}" )
      else()
         set( DBS_OPERATING_SYSTEM "${win_prod_name} ${win_sp}" )
      endif()
      
      execute_process(
         COMMAND ${CMAKE_CXX_COMPILER}
         ERROR_VARIABLE tmp
         OUTPUT_STRIP_TRAILING_WHITESPACE
         OUTPUT_QUIET
         )
      string( REGEX REPLACE "^(.*)Copyright.*" "\\1"
         tmp "${tmp}" )
      string( STRIP "${tmp}" DBS_CXX_COMPILER_VER )
      set( DBS_C_COMPILER_VER "${DBS_CXX_COMPILER_VER}" )
      
      # Did we build 32-bit or 64-bit code?
      execute_process(
         COMMAND ${CMAKE_CXX_COMPILER}
         ERROR_VARIABLE tmp
         OUTPUT_STRIP_TRAILING_WHITESPACE
         OUTPUT_QUIET
         )
      string( REGEX REPLACE ".*for ([0-9x]+)" "\\1"
         tmp "${DBS_CXX_COMPILER_VER}" )
      if( ${tmp} MATCHES "80x86" )
         set( DBS_ISA_MODE "32-bit" )
      elseif( ${tmp} MATCHES "x64" )
         set( DBS_ISA_MODE "64-bit" )
      endif( ${tmp} MATCHES "80x86" )

      # Try to get the build machine name
      execute_process(
         COMMAND "c:/windows/system32/ipconfig.exe" "/all"
         OUTPUT_VARIABLE windows_ip_configuration
         )
      string( REGEX REPLACE ".*Host Name[.: ]+([A-z]+).*Primary.*" "\\1" DBS_TARGET "${windows_ip_configuration}" )
      
      # end elseif( WIN32 )
   endif()

endmacro()

#------------------------------------------------------------------------------#
# End 
#------------------------------------------------------------------------------#
