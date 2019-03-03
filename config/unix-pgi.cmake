#-----------------------------*-cmake-*----------------------------------------#
# file   config/unix-pgi.cmake
# brief  Establish flags for Linux64 - Portland Group C/C++
# note   Copyright (C) 2016-2019 Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

# History
# ----------------------------------------
# 6/13/2016  - IPO settings moved to compilerEnv.cmake
#              (CMAKE_INTERPROCEDURAL_OPTIMIZATION=ON).

set( draco_isPGI 1 CACHE BOOL "Are we using PGI CXX? -> ds++/config.h" )

# NOTE: You may need to set TMPDIR to a location that the compiler can
#       write temporary files to.  Sometimes the default space on HPC
#       is not sufficient.

#
# Sanity Checks
#

if( BUILD_SHARED_LIBS )
  message( FATAL_ERROR "Feature not available - yell at KT." )
endif( BUILD_SHARED_LIBS )

if( CMAKE_CXX_COMPILER_VERSION VERSION_LESS 13.0 )
  message( FATAL_ERROR "PGI must be version 13.0 or later to support C++ features used by this software.")
endif()

#
# Compiler Flags
#

# [2015-03-09 KT] Random123 has a lot of code that produces the warning:
#    "warning #111-D:  statement is unreachable."
# Choose to supress these with --diag_suppress 111

if( NOT CXX_FLAGS_INITIALIZED )
  set( CXX_FLAGS_INITIALIZED "yes" CACHE INTERNAL "using draco settings." )

  set( CMAKE_C_FLAGS                "-Kieee -Mdaz -pgf90libs" )
  set( CMAKE_C_FLAGS_DEBUG          "-g -O0")
  set( CMAKE_C_FLAGS_RELEASE        "-O3 -DNDEBUG" )
  set( CMAKE_C_FLAGS_MINSIZEREL     "${CMAKE_C_FLAGS_RELEASE}" )
  set( CMAKE_C_FLAGS_RELWITHDEBINFO "-O3 -gopt" )
  set( CMAKE_CXX_FLAGS              "${CMAKE_C_FLAGS} --c++11 --no_using_std --no_implicit_include --display_error_number --diag_suppress 940 --diag_suppress 11 --diag_suppress 450 --diag_suppress 111 -DNO_PGI_OFFSET" )

  # Extra flags for pgCC-11.2+
  # --nozc_eh    (default for 11.2+) Use low cost exception handling. This
  #              option appears to break our exception handling model resulting
  #              in SEGV.
  # This may be related to PGI bug 1858 (http://www.pgroup.com/support/release_tprs.htm).

  set( CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_C_FLAGS_DEBUG}")
  set( CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_C_FLAGS_RELEASE} -Munroll=c:10 -Mautoinline=levels:10 -Mvect=sse -Mflushz -Mlre")
  set( CMAKE_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_RELEASE}")
  set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -gopt" )

  # Use the C99 standard
  set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -c99")

endif()

##---------------------------------------------------------------------------##
# Ensure cache values always match current selection
##---------------------------------------------------------------------------##
set( CMAKE_C_FLAGS                "${CMAKE_C_FLAGS}"                CACHE STRING "compiler flags" FORCE )
set( CMAKE_C_FLAGS_DEBUG          "${CMAKE_C_FLAGS_DEBUG}"          CACHE STRING "compiler flags" FORCE )
set( CMAKE_C_FLAGS_RELEASE        "${CMAKE_C_FLAGS_RELEASE}"        CACHE STRING "compiler flags" FORCE )
set( CMAKE_C_FLAGS_MINSIZEREL     "${CMAKE_C_FLAGS_MINSIZEREL}"     CACHE STRING "compiler flags" FORCE )
set( CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}" CACHE STRING "compiler flags" FORCE )

set( CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS}"                CACHE STRING "compiler flags" FORCE )
set( CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG}"          CACHE STRING "compiler flags" FORCE )
set( CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE}"        CACHE STRING "compiler flags" FORCE )
set( CMAKE_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_MINSIZEREL}"     CACHE STRING "compiler flags" FORCE )
set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}" CACHE STRING "compiler flags" FORCE )

toggle_compiler_flag( OPENMP_FOUND ${OpenMP_C_FLAGS} "C;CXX;EXE_LINKER" "" )

#------------------------------------------------------------------------------#
# End config/unix-pgi.cmake
#------------------------------------------------------------------------------#
