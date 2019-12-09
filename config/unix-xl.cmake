#-----------------------------*-cmake-*----------------------------------------#
# file   config/unix-xl.cmake
# author Gabriel Rockefeller
# date   2012 Nov 1
# brief  Establish flags for Linux64 - IBM XL C++
# note   Copyright (C) 2016-2019 Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

# Ref:
# https://www.ibm.com/support/knowledgecenter/SSXVZZ_16.1.0/com.ibm.xlcpp161.lelinux.doc/compiler_ref/opt_langlvl.html

#
# Compiler flag checks
#
include(platform_checks)
query_openmp_availability()

#
# Compiler Flags
#

if( NOT CXX_FLAGS_INITIALIZED )
   set( CXX_FLAGS_INITIALIZED "yes" CACHE INTERNAL "using draco settings." )

  # On Darwin, we also need this config file:
  # -F/projects/opt/ppc64le/ibm/xlc-16.1.1.2/xlC/16.1.1/etc/xlc.cfg.rhel.7.5.gcc.7.3.0.cuda.9.2
  # -qfloat=nomaf -qxlcompatmacros
  set( CMAKE_C_FLAGS                "-g" ) # -qarch=auto
  # Sequoia
  if( CMAKE_CXX_COMPILER_VERSION VERSION_LESS 13.0 )
    string( APPEND CMAKE_C_FLAGS " -qinfo=all -qflags=i:w -qsuppress=1540-0072")
    string( APPEND CMAKE_C_FLAGS " -qsuppress=1506-1197" )
  endif()
  # 2019-04-03 IBM support asks that we not use '-qcheck' due to compiler issues.
  set( CMAKE_C_FLAGS_DEBUG          "-O0 -qsmp=omp:noopt -qoffload -qfullpath -DDEBUG") # -qcheck
  set( CMAKE_C_FLAGS_RELWITHDEBINFO
    "-O3 -qhot=novector -qsimd=auto -qstrict=nans:operationprecision" )
  set( CMAKE_C_FLAGS_RELEASE        "${CMAKE_C_FLAGS_RELWITHDEBINFO} -DNDEBUG" )
  set( CMAKE_C_FLAGS_MINSIZEREL     "${CMAKE_C_FLAGS_RELEASE}" )

   # Email from Roy Musselman <roymuss@us.ibm.com, 2019-03-21:
   # For C++14, add -qxflag=disable__cplusplusOverride
   set( CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -qxflag=disable__cplusplusOverride -Wno-undefined-var-template")
   set( CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_C_FLAGS_DEBUG}")
   set( CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_C_FLAGS_RELEASE}")
   set( CMAKE_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_RELEASE}")
   set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}" )

  if( ${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "ppc64le")
    string( APPEND CMAKE_C_FLAGS " -qarch=pwr9 -qtune=pwr9" )
  endif()

endif()

##---------------------------------------------------------------------------##
# Ensure cache values always match current selection
##---------------------------------------------------------------------------##
set( CMAKE_C_FLAGS                "${CMAKE_C_FLAGS}"                CACHE STRING
  "compiler flags" FORCE )
set( CMAKE_C_FLAGS_DEBUG          "${CMAKE_C_FLAGS_DEBUG}"          CACHE STRING
  "compiler flags" FORCE )
set( CMAKE_C_FLAGS_RELEASE        "${CMAKE_C_FLAGS_RELEASE}"        CACHE STRING
  "compiler flags" FORCE )
set( CMAKE_C_FLAGS_MINSIZEREL     "${CMAKE_C_FLAGS_MINSIZEREL}"     CACHE STRING
  "compiler flags" FORCE )
set( CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}" CACHE STRING
  "compiler flags" FORCE )

set( CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS}"                CACHE
  STRING "compiler flags" FORCE )
set( CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG}"          CACHE
  STRING "compiler flags" FORCE )
set( CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE}"        CACHE
  STRING "compiler flags" FORCE )
set( CMAKE_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_MINSIZEREL}"     CACHE
  STRING "compiler flags" FORCE )
set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}" CACHE
  STRING "compiler flags" FORCE )

#toggle_compiler_flag( DRACO_SHARED_LIBS "-qnostaticlink" "EXE_LINKER" "")

# CMake will set OpenMP_C_FLAGS to '-qsmp.'  This option turns on
# OpenMP but also activates the auto-parallelizer.  We don't want to
# enable the 2nd feature so we need to specify the OpenMP flag to be
# '-qsmp=omp.'
if( CMAKE_CXX_COMPILER_VERSION VERSION_LESS 13.0 )
  toggle_compiler_flag( OPENMP_FOUND             "-qsmp=omp" "C;CXX;EXE_LINKER"
    "" )
#else()
  # toggle_compiler_flag( OPENMP_FOUND             "-qsmp=noauto"
  # "C;CXX;EXE_LINKER" "" )
endif()

#------------------------------------------------------------------------------#
# End config/unix-xl.cmake
#------------------------------------------------------------------------------#
